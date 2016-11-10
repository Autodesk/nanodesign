#!/usr/bin/env python
""" 
This module is used to create a DNA structure from a caDNAno DNA origami design file. 

The virtual helices from the CadnanoDesign object are used to create a list of DnaStructureHelix objects. 
DnaBase objects are created from the scaffold and staple bases defined for each virtual helix and stored 
in the appropriate DnaStructureHelix object. A DnaStructure object object is created for the design. It 
stores the list of DnaStructureHelix objects and a list of all the bases defined for the design. 

This code is initially based on a direct translation of the set MATLAB scripts to convert a caDNAno design 
to a CanDo .cndo file from Mark Bathe's Laboratory for Computational Biology & Biophysics at MIT.
"""
from collections import OrderedDict
import csv
import logging
import os
import re
import string
import sys
import time
import numpy as np
from math import sqrt
from ..dna_sequence_data import dna_sequence_data

from .design import CadnanoDesign,CadnanoVirtualHelix,CadnanoBase
from .reader import CadnanoReader 
from .common import CadnanoLatticeType
from .utils import generate_coordinates,get_start_coordinates_angle,vrrotvec2mat,deg2rad,bp_interp,find_row

from ...data.base import DnaBase 
from ...data.strand import DnaStrand
from ...data.dna_structure import DnaStructure,DnaStructureHelix
from ...data.lattice import Lattice,SquareLattice,HoneycombLattice
from ...data.parameters import DnaPolarity,DnaParameters

class StrandType:
    SCAFFOLD = 0
    STAPLE   = 1

class CadnanoConvertDesign(object):
    """ The CadnanoConvertDesign class is used to create a DNA structure from a caDNAno DNA origami design. 

        Attributes:
            base_id (int): The current base ID used when creating a DnaBase object.
            base_map (OrderedDict[Tuple]): The dict that stores DnaBase objects using the tuple (helix_num,helix_pos) as a key.
            dna_parameters (DnaParameters): The DNA parameters to use when creating the 3D geometry for the design.
            dna_structure (DnaStructure): The DnaStructure object. 
            staple_colors (List[StapleColor]: The list of caDNAno staple strand colors.
    """
    def __init__(self, dna_parameters):
        """ Initialize the CadnanoConvertDesign object.

            Arguments:
                dna_parameters (DnaParameters): The DNA parameters to use when creating the 3D geometry for the design.
        """
        self._logging_level = logging.INFO
        self._setup_logging()
        self.dna_structure = None 
        self.staple_colors = []
        self.dna_parameters = dna_parameters
        self.base_id = 0
        self.base_map = OrderedDict()

    def _set_logging_level(self,level):
        """Set logging level."""
        self._logger.setLevel(level)

    def _setup_logging(self):
        """ Set up logging. """
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(self._logging_level)
        # Create a console handler and set format.
        if not len(self._logger.handlers):
            console_handler = logging.StreamHandler()
            formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
            console_handler.setFormatter(formatter)
            self._logger.addHandler(console_handler)

    def _add_staple_color(self, staple_color, vhelix_num):
        """ Add a caDNAno staple color for the given virtual helix. 

            Arguments:
                staple_color (List[int]): The caDNAno color definition consisting of a list of two integer
                    values, position and color value.
                vhelix_num (int): The caDNAno virtual helix number. 
        """
        self.staple_colors.append(self.StapleColor(staple_color,vhelix_num))

    def _get_base_index(self, num, pos, strand_type):
        """ Create a tuple (helix_num,helix_pos,strand_type) used as a key to store bases in the base_map dict.

            Arguments:
                num (int): The base virtual helix number. 
                pos (int): The base virtual helix position. 
                strand_type (StrandType): The type of strand scaffold or staple the base is part of.

            Returns The tuple (helix_num,helix_pos,strand_type) used to uniquely identify a base.
         """
        base_index = (int(num), int(pos), int(strand_type))
        return base_index 

    def _get_base(self, base_index):
        """ Get a base from the base_map dict. If the base is not in the dict then create and add it.

            Arguments:
                base_index (Tuple): The tuple (helix_num,helix_pos,strand_type) used to uniquely identify a base.

            Returns a DnaBase object.
        """
        if base_index not in self.base_map:
            base = DnaBase(self.base_id)
            self.base_map[base_index] = base
            self.base_id += 1
        else:
            base = self.base_map[base_index]
        return base
    #__def _get_base

    def create_structure(self, design, modify=False):
        """ Create a DNA structure from a caDNAno DNA origami design.

            Arguments:
                design (CadnanoDesign): A caDNAno DNA origami design.
                modify (bool): If true then modify the caDNAno DNA origami design for inserts and deletes.

            Returns:
                DnaStructure: The DNA structure object containing topological and geometric information for the
                              caDNAno DNA origami design model.

            The virtual helices from the CadnanoDesign object are design are used to create a list of 
            DnaStructureHelix objects. DnaBase objects are created from the scaffold and staple bases 
            defined for each virtual helix and stored in the appropriate DnaStructureHelix object. The
            global list of DnaBase objects are stored in the self.base_map dict.
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.info("Distance between adjacent helices %g " % self.dna_parameters.helix_distance)
        self._logger.info("Helix radius %g " % self.dna_parameters.helix_radius) 
        # Reset these in case this function is called multiple times. 
        self.base_id = 0
        self.base_map = OrderedDict()

        # Create a list of DnaStructureHelix objects for the design. 
        helices = self._create_structure_topology_and_geometry(design)
        self._logger.info("Number of bases in design %d " % len(self.base_map))

        # Set the bases up and down attributes pointing to terminal bases to point to None. 
        # These scaffold and staple terminal bases are automatically added when adding caDNAno 
        # virtual helix bases.
        print_dna_model_base_map = False
        num_bases = 0
        n = 0
        term_bases = []
        if (-1,-1,1) in self.base_map:
            term_bases.append(self.base_map[(-1,-1,1)].id)
        if (-1,-1,0) in self.base_map:
            term_bases.append(self.base_map[(-1,-1,0)].id)
        for bindex, base in self.base_map.items():
            if base.id not in term_bases:
                num_bases += 1
            if base.up and base.up.id in term_bases:
               base.up = None
            if base.down and base.down.id in term_bases:
               base.down = None
            up = base.up.id if base.up else -1
            down = base.down.id if base.down else -1
            across = base.across.id if base.across else -1
            if print_dna_model_base_map:
                self._logger("%4d: bindex %12s  base id %4d  up %4s  down %4s  across %4s  h %4d  p %4d  scaf %d" %
                    ( n, bindex, base.id, up, down, across, base.h, base.p, base.is_scaf))
            n += 1
        #__for bindex, base in dna_model.base_map.items():

        # Create an array of all of the DnaBase objects created for the design.
        base_connectivity = [None]*num_bases
        num_bases = 0
        for bindex, base in self.base_map.items():
            if base.id not in term_bases:
                base.id = num_bases
                base_connectivity[base.id] = base
                num_bases += 1
        #_for bindex, base in self.base_map.items()
        if self._logger.getEffectiveLevel() == logging.DEBUG:
            self._logger.debug("---------- base_connectivity  ---------- ")
            self._logger.debug("size of base_connectivity: %d " % len(base_connectivity))
            for i in xrange(0,num_bases):
                base = base_connectivity[i]
                up = base.up.id if base.up else -1
                down = base.down.id if base.down else -1
                across = base.across.id if base.across else -1
                self._logger.debug("%4d  id %4d  h %4d  p %4d  up %4d  down %4d  across %4d  scaf %d" %
                    ( i, base.id, base.h, base.p, up, down, across, base.is_scaf))
        #__if print_base_connectivity_p
        #__if self._logger.getEffectiveLevel() == logging.DEBUG

        # Remove deleted bases.
        if (modify):
            self._delete_bases(helices, base_connectivity)
            print_base_connectivity = False
            if print_base_connectivity:
                self._logger.info("Size of topology after deletes %d" % len(base_connectivity))
                self._logger.info("---------- base_connectivity after deletes ---------- ")
                self._logger.info("size of base_connectivity: %d " % len(base_connectivity))
                for i,base in enumerate(base_connectivity): 
                    up = base.up.id if base.up else -1
                    down = base.down.id if base.down else -1
                    across = base.across.id if base.across else -1
                    self._logger.info("%4d: base id %4d  up %4s  down %4s  across %4s  h %4d  p %4d  scaf %d" %
                        ( i, base.id, up, down, across, base.h, base.p, base.is_scaf))
            #__if print_base_connectivity_p
        #__if (modify)

        # Add inserted bases.
        if (modify):
            self._insert_bases(helices, base_connectivity)
            print_base_connectivity = False
            if print_base_connectivity:
                self._logger.info("Size of topology after inserts %d" % len(base_connectivity))
                self._logger.info("---------- base_connectivity after inserts ---------- ")
                for i,base in enumerate(base_connectivity):
                    up = base.up.id if base.up else -1
                    down = base.down.id if base.down else -1
                    across = base.across.id if base.across else -1
                    self._logger.info("%4d: base id %4d  up %4s  down %4s  across %4s  h %4d  p %4d  scaf %d" %
                        ( i, base.id, up, down, across, base.h, base.p, base.is_scaf))
            #__if print_base_connectivity_p
        #__if (modify)

        # Create a DnaStructure object to store the base connectivty and helices.
        name = "dna structure"
        self.dna_structure = DnaStructure(name, base_connectivity, helices, self.dna_parameters)
        self.dna_structure.set_lattice_type(design.lattice_type)

        # Generate strands.
        self.dna_structure.create_strands()
        strands = self.dna_structure.strands
        if strands == None:
            self._logger.error("Create strands failed.")
            sys.exit(1)
        self._set_strands_colors(strands)
        self.dna_structure.strands = strands
        self._logger.info("Number of strands %d " % len(strands)) 
        if self._logger.getEffectiveLevel() == logging.DEBUG:
            for strand in strands:
                self._logger.debug("Strand %4d  bases %4d  scaf %6s  start helix %4d  pos %4d" % (strand.id, 
                    len(strand.tour), strand.is_scaffold, strand.tour[0].h,  strand.tour[0].p))
        #__if self._logger.getEffectiveLevel() == logging.DEBUG

        # Calculate staple ends.
        self.dna_structure.staple_ends = self._calculate_staple_ends(strands)

        # Set possible cross-overs. 
        self._set_possible_crossovers(design)
 
        return self.dna_structure
    #__def create_structure

    def _delete_bases(self, helices, base_connectivity):
        """ Remove bases from helices and base_connectivity.  """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("==================== delete bases ====================")

        # Delete bases from each helix.
        num_deleted_bases = 0
        for helix in helices:
            deleted_bases = helix.process_base_deletes()
            num_deleted_bases += len(deleted_bases)
            # Remove bases from base_connectivity.
            for base in deleted_bases:
                self._logger.debug("Delete base %d" % base.id)
                base_connectivity.remove(base)
            #__for base in deleted_bases
        #__for helix in helices

        # Renumber base IDs.
        if num_deleted_bases != 0:
            self._logger.info("Number of deleted bases %d" % num_deleted_bases) 
            self._renumber_baseIDs(base_connectivity)

    def _insert_bases(self, helices, base_connectivity):
        """ Insert bases into helices and base_connectivity[]. 

            Arguments:
                helices(List[DnaStructureHelix]): The list of helices for the structure.
                base_connectivity (List[DnaBase]): The list of DNA bases for the structure.

           The DnaBase.num_inserts attribute determines the number of bases inserted at that base. 
           Bases are inserted in the 3' direction. 
        """
        self._logger.setLevel(logging.INFO)
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("==================== insert bases ====================")
        num_bases = len(base_connectivity)
        base_id = num_bases
        y_up_vec   = np.array([0,  1.0, 0], dtype=float)
        y_down_vec = np.array([0, -1.0, 0], dtype=float)

        # Create a list of inserts from helices base lists.
        num_insert = 0
        num_insert_bases = 0
        base_inserts = OrderedDict()
        helix_map = {}  # Stores DnaStructureHelix objects used later to update base lists.
        for helix in helices:
            helix_map[helix.id] = helix
            for base in helix.staple_bases:
                if base.num_insertions != 0:
                    num_insert_bases += base.num_insertions 
                    base_inserts[base.id] = base
            #__for base in helix.staple_bases
            for base in helix.scaffold_bases:
                if base.num_insertions != 0:
                    base_inserts[base.id] = base
                    self._logger.debug("Insert base id %d  h %d  pos %d" % (base.id, base.h, base.p))
            #__for base in helix.scaffold_bases
        #__for helix in helices
        self._logger.debug("Number of base inserts %d" % len(base_inserts))
        if len(base_inserts) == 0:
            return 

        # Iterated over bases with inserts. 
        processed_bases = set()
        new_bases = {} # Stores lists of new insert bases by helix ID.
        for curr_base in base_inserts.values():
            if curr_base.id in processed_bases:
                continue
            processed_bases.add(curr_base.id)
            coords = curr_base.coordinates
            self._logger.debug("-------------------------------------------")
            self._logger.debug("Current base      ID %3d  h %3d  p %3d  (%6.2f %6.2f %6.2f)" % 
                (curr_base.id, curr_base.h, curr_base.p, coords[0], coords[1], coords[2]))
            if curr_base.up:
                up_base = curr_base.up
                coords = up_base.coordinates
                self._logger.debug("Current base up   ID %3d  h %3d  p %3d  (%6.2f %6.2f %6.2f) " % 
                    (up_base.id, up_base.h, up_base.p, coords[0], coords[1], coords[2]))
            if curr_base.down:
                down_base = curr_base.down
                coords = down_base.coordinates
                self._logger.debug("Current base down ID %3d  h %3d  p %3d  (%6.2f %6.2f %6.2f) " % 
                    (down_base.id, down_base.h, down_base.p, coords[0], coords[1], coords[2]))
            helix_axis = curr_base.ref_frame[:,2]
            curr_across = curr_base.across

            if curr_across != None:
                self._logger.debug("Current across base ID %d" % curr_across.id)
                neighbor_across_up = curr_across.up
                processed_bases.add(curr_across.id)
                # Make sure that the 3'-neighbor of the current base is to the right.
                tmp_a = curr_base.ref_frame[:,2]
                if (np.linalg.norm(tmp_a - y_up_vec) < 1e-10):  # Reference axis e3 poins to the right.
                    self._logger.debug("Reference axis points to right")
                    # If curr_base is not a preferred base then swap with its across base.
                    if not curr_base.is_scaf: 
                        tmp_base = curr_base
                        curr_base  = curr_across
                        curr_across = tmp_base
                    #__if not curr_base.is_scaf
                elif (np.linalg.norm(tmp_a - y_down_vec) < 1e-10): # Reference axis e3 points to the left.
                    self._logger.debug("Reference axis points to left ")
                    # If curr_base is a preferred base then swap with its across base.
                    if curr_base.is_scaf: 
                        tmp_base = curr_base
                        curr_base  = curr_across
                        curr_across = tmp_base
                    #__if curr_base.is_scaf
                else: 
                    self._logger.error("Can't find reference axis when inserting bases.")
                    sys.exit(0)
                #__if (np.linalg.norm(tmp_a - y_up_vec) < 1e-10) 
            #__if curr_across != None

            # Insert base into dsDNA. 
            if curr_across != None:
                base_id = self._insert_bases_dsDNA(curr_base, base_id, helix_axis, new_bases, base_connectivity)

            # Insert base into ssDNA. 
            else:
                base_id = self._insert_bases_ssDNA(curr_base, base_id, helix_axis, new_bases, base_connectivity)
            #__if curr_across != None
        #__for base in inserted_bases

        # Add the inserted bases into their parent helices.
        for helix_id in new_bases.keys():
            base_list = new_bases[helix_id]
            self._logger.debug("Helix ID %d  add num bases %d " % (helix_id, len(base_list)))
            helix_map[helix_id].insert_bases(base_list)
        #__for helix_id in new_bases.keys()

        self._logger.info("Number of inserted bases %d" % num_insert_bases)
    #__def _insert_bases

    def _insert_bases_ssDNA(self, curr_base, base_id, helix_axis, new_bases, base_connectivity):
        """ Insert a number of bases into ssDNA. 

            Arguments:
                curr_base (DnaBase): The base where insertaions are made.
                base_id (int): The current base ID. This is incremented to set new insert base IDs.
                helix_axis (NumPy 3 ndarray[float]): The helix axis.
                new_bases (Dict): The map of helix IDs to lists of new insert bases.
                base_connectivity (List[DnaBase]): The list of DNA bases for the structure.

            Returns base_id. 

            New bases are inserted in the 3' direction for curr_base. The new DnaBase objects are stored
            in new_bases{} by helix ID and appended to base_connectivity[].
        """
        dist_bp = self.dna_parameters.base_pair_rise         # Rise between two neighboring base-pairs (nm).
        ang_bp = self.dna_parameters.base_pair_twist_angle   # Twisting angle between two neighboring base-pairs (degree)
        dy = np.array([0, 0.1, 0], dtype=float)
        y_up_vec   = np.array([0,  1.0, 0], dtype=float)
        y_down_vec = np.array([0, -1.0, 0], dtype=float)

        # Find the neighboring bases.
        curr_across = curr_base.across
        neighbor_down = curr_base.down
        if curr_across != None:
            neighbor_across_up = curr_across.up
        else:
            neighbor_across_up = None

        # Interpolate positions and orientations.
        num_inserts = curr_base.num_insertions
        curr_coords = curr_base.coordinates
        curr_frame = curr_base.ref_frame.copy()
        curr_frame = curr_frame.reshape(3,3)
        next_coords = curr_coords + [0, -dist_bp*helix_axis[1], 0]
        rot_mat = vrrotvec2mat(y_up_vec, deg2rad(ang_bp))
        next_frame = np.dot(rot_mat,curr_frame)
        [insert_coords, insert_frames] = bp_interp(curr_coords, curr_frame, next_coords, next_frame, num_inserts)
        last_base = None
        self._logger.debug("Current coords %g %g %g" % (curr_coords[0], curr_coords[1], curr_coords[2]))
        self._logger.debug("Next coords %g %g %g" % (next_coords[0], next_coords[1], next_coords[2]))
        self._logger.debug("Insert ssDNA")
        self._logger.debug("num_inserts %d " % num_inserts)

        # Create bases to insert.
        for i in xrange(0,num_inserts):
            base = DnaBase(base_id)
            if curr_base.h not in new_bases:
                new_bases[curr_base.h] = []
            new_bases[curr_base.h].append(base)
            base_id += 1
            base_connectivity.append(base)

            # Set base connectivity.
            if (i == 0):
                curr_base.down = base 
                base.up = curr_base
            else:
                base.up = last_base 
                last_base.down = base

            if (i == num_inserts-1):
                if neighbor_down != None:
                    neighbor_down.up = base
                base.down = neighbor_down
            else:
                base.down = last_base
                #last_base.up = base

            self._logger.debug("Insert %d coords %g %g %g" % (i, insert_coords[i,0], insert_coords[i,1], insert_coords[i,2]))
            base.h = curr_base.h
            base.p = curr_base.p
            base.is_scaf = curr_base.is_scaf
            base.nt_coords = curr_base.nt_coords + dy*(i+1+1)/2
            base.coordinates = insert_coords[i]
            base.ref_frame = insert_frames[:,:,i]
            last_base = base
        #__for i in xrange(0,num_inserts)

        return base_id 
    #__def _insert_bases_ssDNA

    def _insert_bases_dsDNA(self, curr_base, base_id, helix_axis, new_bases, base_connectivity):
        """ Insert a number of bases into dsDNA. 

            Arguments:
                curr_base (DnaBase): The base where insertaions are made. 
                base_id (int): The current base ID. This is incremented to set new insert base IDs. 
                helix_axis (NumPy 3 ndarray[float]): The helix axis. 
                new_bases (Dict): The map of helix IDs to lists of new insert bases. 
                base_connectivity (List[DnaBase]): The list of DNA bases for the structure.

            Returns base_id. 

            New bases are inserted in the 3' direction for curr_base. The new DnaBase objects are stored
            in new_bases{} by helix ID and appended to base_connectivity[].
        """
        dist_bp = self.dna_parameters.base_pair_rise         # Rise between two neighboring base-pairs (nm).
        ang_bp = self.dna_parameters.base_pair_twist_angle   # Twisting angle between two neighboring base-pairs (degree)
        dy = np.array([0, 0.1, 0], dtype=float)
        y_up_vec   = np.array([0,  1.0, 0], dtype=float)
        y_down_vec = np.array([0, -1.0, 0], dtype=float)

        # Find the neighboring bases.
        neighbor_down = curr_base.down
        curr_across = curr_base.across
        if curr_across != None:
            neighbor_across_up = curr_across.up
        else:
            neighbor_across_up = None

        # Interpolate positions and orientations.
        num_inserts = 2*curr_base.num_insertions
        curr_coords = curr_base.coordinates
        curr_frame = curr_base.ref_frame.copy()
        curr_frame = curr_frame.reshape(3,3)
        #next_coords = curr_coords + [0, dist_bp, 0]
        next_coords = curr_coords + [0, -dist_bp*helix_axis[1], 0]
        rot_mat = vrrotvec2mat(y_up_vec, deg2rad(ang_bp))
        next_frame = np.dot(rot_mat,curr_frame)
        [insert_coords, insert_frames] = _bp_interp(curr_coords, curr_frame, next_coords, next_frame, num_inserts/2)
        self._logger.debug("Insert dsDNA")
        self._logger.debug("num_inserts %d " % num_inserts)
        self._logger.debug("next_coords %s" % str(next_coords))
        self._logger.debug("curr_coords %s" % str(curr_coords))
        self._logger.debug("insert coords %s" % str(insert_coords )) 

        # Create bases to insert.
        last_base1 = None 
        last_base2 = None 
        for i in xrange(0,num_inserts,2):
            base1 = DnaBase(base_id)
            if curr_base.h not in new_bases:
                new_bases[curr_base.h] = []
            new_bases[curr_base.h].append(base1)
            base_id += 1
            base2 = DnaBase(base_id)
            new_bases[curr_base.h].append(base2)
            base_id += 1
            self._logger.debug("Add two bases  Base1 id %d   Base 2 id %d " % (base1.id, base2.id))
            base_connectivity.append(base1)
            base_connectivity.append(base2)

            # Set base connectivity.
            if i == 0:
                curr_base.down = base1
                curr_across.up = base2 
                base1.up = curr_base
                base2.down = curr_across
            else:
                base1.up = last_base1 
                base2.down = last_base2
                last_base1.down = base1
                last_base2.up = base2

            if i == num_inserts-2:
                if neighbor_down != None:
                    neighbor_down.up = base1 
                if neighbor_across_up != None:
                    neighbor_across_up.down = base2 
                base1.down = neighbor_down
                base2.up = neighbor_across_up
            else:
                base1.down = last_base1
                base2.up = last_base2 

            base1.across = base2
            base1.h = curr_base.h
            base1.p = curr_base.p
            base1.is_scaf = curr_base.is_scaf
            base1.nt_coords = curr_base.nt_coords + dy*(i+1+1)/2
            base1.coordinates = insert_coords[i/2]
            base1.ref_frame = insert_frames[i/2]

            base2.across = base1
            base2.h = curr_across.h
            base2.p = curr_across.p
            base2.is_scaf = curr_across.is_scaf
            base2.nt_coords = curr_across.nt_coords + dy*(i+1)/2
            base2.coordinates = insert_coords[i/2]
            base2.ref_frame = insert_frames[i/2]

            last_base1 = base1
            last_base2 = base2
        #__for i in xrange(0,num_inserts,2)

        return base_id 
    #__def _insert_bases_dsDNA

    def _renumber_baseIDs(self, base_connectivity):
        """ Renumber base IDs to be between 0 and num_bases. """
        for i,base in enumerate(base_connectivity): 
            base.id = i

    def _create_structure_topology_and_geometry(self, design):
        """ Create topological and geometrical information for a design. """
        lattice_type = design.lattice_type
        max_vhelix_size = design.max_base_id+1
        row_list = []
        col_list = []
        structure_helices = [] 
        vhelices = design.helices
        self._logger.setLevel(logging.INFO)
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("==================== create structure topology and geometry ====================")

        for i,vhelix in enumerate(vhelices):
            self._logger.debug("---------- process virtual helix num %d ----------" % vhelix.num );
            row = vhelix.row 
            col = vhelix.col 
            num = vhelix.num 
            row_list.append(row)
            col_list.append(col)
            self._logger.debug("num: %d " % num)
            self._logger.debug("row: %d " % row)
            self._logger.debug("col: %d " % col)

            if ( num % 2 == 0 ):
                scaffold_polarity = DnaPolarity.FIVE_PRIME
                self._logger.debug("scaffold polarity 5' to 3'")
            else:
                scaffold_polarity = DnaPolarity.THREE_PRIME
                self._logger.debug("scaffold polarity 3' to 5'")

            # Set staple colors
            stap_colors = vhelix.staple_colors 
            for color in stap_colors:
                self._add_staple_color(color, num)

            # Get the bases for the helix.
            scaffold_bases, staple_bases = self._create_single_helix(vhelix, lattice_type)
            self._logger.debug("Number of scaffold bases %d " % len(scaffold_bases)) 
	    self._logger.debug("Number of staple bases %d " % len(staple_bases)) 
            if self._logger.getEffectiveLevel() == logging.DEBUG:
                s = ""
                for base in scaffold_bases:
                    s += str(base.p) + " "
                self._logger.debug("Scaffold bases positions %s " % s) 
                s = ""
                for base in staple_bases:
                    s += str(base.p) + " "
                self._logger.debug("Staple bases positions %s " % s) 
   
            # Generate the helix axis coordinates and frames, and DNA helix nucleotide coordinates.
            axis_coords, axis_frames, scaffold_coords, staple_coords = \
                generate_coordinates(self.dna_parameters, lattice_type, row, col, num, scaffold_bases, staple_bases)

            # Create a dna structure object that stores the helix information. 
            structure_helix = DnaStructureHelix(i, num, scaffold_polarity, axis_coords, axis_frames, 
                scaffold_coords, staple_coords, scaffold_bases, staple_bases)
            structure_helix.lattice_num = num
            structure_helix.lattice_row = row
            structure_helix.lattice_col = col
            structure_helix.lattice_max_vhelix_size = max_vhelix_size 
            structure_helices.append(structure_helix)
        #__for vhelix in vhelices
        return structure_helices
    #__def _create_structure_topology_and_geometry

    def _create_single_helix(self, vhelix, lattice_type):
        """ Create the geometry for a single cadnano virtual helix.

            Arguments:
                vhelix (CadnanoVirtualHelix): A cadnano virtual helix.
                lattice_type (CadnanoLatticeType): Cadnano lattice type. 

            Returns: 
                scaffold_bases (List[DnaBase]): The list of scaffold bases defined for the helix. 
                staple_bases (List[DnaBase]): The list of staple bases defined for the helix. 
        """
        #self._logger.debug("------------------- _create_single_helix ------------------- " )
        scaffolds = vhelix.scaffold_strands 
        staples = vhelix.staple_strands 
        deletions = vhelix.deletions 
        insertions = vhelix.insertions 
        row = vhelix.row 
        col = vhelix.col 
        num = vhelix.num 
        num_bases = len(scaffolds)

        # Iterate over the vhelix scaffold and staple bases.
        scaffold_bases = []
        staple_bases = []
        for i in xrange(0,num_bases):
            current_scaffold = scaffolds[i] 
            current_staple = staples[i] 
            helix_pos = i

            # If the base exists in the scaffold strand.
            if (current_scaffold.initial_strand >= 0) or (current_scaffold.final_strand >= 0):
                base = self._add_base(StrandType.SCAFFOLD, current_scaffold, StrandType.STAPLE, current_staple, helix_pos, num)
                base.num_deletions = deletions[i]
                base.num_insertions = insertions[i]
                scaffold_bases.append(base)

            # if the base exists in the staple strand
            if ((current_staple.initial_strand >= 0) or (current_staple.final_strand >= 0)):
                base = self._add_base(StrandType.STAPLE, current_staple, StrandType.SCAFFOLD, current_scaffold, helix_pos, num)
                staple_bases.append(base)
                base.num_deletions = deletions[i]
                base.num_insertions = insertions[i]

        #__for i in xrange(0,num_lattice)__

        return scaffold_bases, staple_bases 

    def _add_base(self, base_type, base, paired_base_type, paired_base, helix_pos, helix_num):
        """ Create a base from a cadnano base for a given location in a virtual helix. 
 
            Arguments:
                base_type (StrandType): The type of helix strand, SCAFFOLD or STAPLE, the cadnano base is part of. 
                base (CadnanoBase): The cadnano base. 
                paired_base_type (StrandType): The type of helix paired strand, SCAFFOLD or STAPLE, the cadnano 
                    paired base is part of. 
                paired_base (CadnanoBase): The cadnano paired base. 
                helix_pos (int): The position of the base in the virtual helix.
                helix_num (int): The number of the virtual helix.
        """
        # Create the helix base. 
        base_index = self._get_base_index(helix_num, helix_pos, base_type) 
        new_base = self._get_base(base_index)
        new_base.is_scaf = (base_type == StrandType.SCAFFOLD)
        new_base.h = helix_num
        new_base.p = helix_pos

        #  Add the 5'-neighbor.
        base_index = self._get_base_index(base.initial_strand, base.initial_base, base_type) 
        five_base = self._get_base(base_index)
        new_base.up = five_base

        # Add the 3'-neighbor
        base_index = self._get_base_index(base.final_strand, base.final_base, base_type) 
        three_base = self._get_base(base_index)
        new_base.down = three_base

        # Watson-Crick neighbor.
        if ( (paired_base.initial_strand >= 0) or (paired_base.final_strand >= 0)):
            base_index = self._get_base_index(helix_num, helix_pos, not base_type)
            wc_base = self._get_base(base_index)
            new_base.across = wc_base
        else:
            new_base.across = None

        return new_base
    #__def _add_base

    def _calculate_staple_ends(self, strands):
        """ Find the start helix and position for the staple strands in the structure. 

            Arguments:
                strands (List[DnaStrand]): The list a DnaStrand objects.

            Returns a Dict mapping 2-tuples (start helix, start position) to strand ID.
        """
        num_strands = len(strands)
        staple_ends = {}

        for i in xrange(0,num_strands):
            strand = strands[i]
            is_scaf_strand = strand.tour[0].is_scaf
            for j in xrange(0,len(strand.tour)):
                is_scaf_base = strand.tour[j].is_scaf
            if (not is_scaf_strand):
                h0 = strand.tour[0].h
                p0 = strand.tour[0].p
                h1 = strand.tour[-1].h
                p1 = strand.tour[-1].p
                staple_ends[(h0,p0)] = i+1
        #__for i in xrange(0,num_strands)__
        return staple_ends
    #__def _calculate_staple_ends

    def _set_possible_crossovers(self,design):
        """ Set the possible cross-overs for scaffold and staple strands.
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.setLevel(logging.INFO)
        self._logger.debug("-------------------- set_possible_crossovers --------------------")
        lattice_type = design.lattice_type
        dist_bp = self.dna_parameters.base_pair_rise 
        scaffold_crossovers = []
        staple_crossovers = [] 
        structure_helices_coord_map = self.dna_structure.structure_helices_coord_map
        for vhelix in design.helices:
            num = vhelix.num 
            col = vhelix.col
            row = vhelix.row
            init_coord,_ = get_start_coordinates_angle(self.dna_parameters, lattice_type, row, col, num)
            self._logger.debug(">>> vhelix: num: %d  row: %d  col: %d " % (num, row, col))
            staple_crossovers = vhelix.possible_staple_crossovers
            scaffold_crossovers = vhelix.possible_scaffold_crossovers
            self._logger.debug("            num staple cross-overs: %d " % len(staple_crossovers)) 
            shelix = structure_helices_coord_map[(row,col)]
            for cross_vh,index in staple_crossovers:
                self._logger.debug("            staple cross-over: %d,%d " % (cross_vh.num,index))
                cross_sh = structure_helices_coord_map[(cross_vh.row,cross_vh.col)]
                coord = init_coord + np.array([0, dist_bp*index, 0]);
                shelix.possible_staple_crossovers.append((cross_sh,index,coord))
            self._logger.debug("            num scaffold cross-overs: %d " % len(scaffold_crossovers)) 
            for cross_vh,index in scaffold_crossovers:
                self._logger.debug("            scaffold cross-over: %d,%d " % (cross_vh.num,index))
                cross_sh = structure_helices_coord_map[(cross_vh.row,cross_vh.col)]
                coord = init_coord + np.array([0, dist_bp*index, 0]);
                shelix.possible_scaffold_crossovers.append((cross_sh,index,coord))
        #__for vhelix in design.helices

    def _set_strands_colors(self, strands):
        """ Set the color for staple strands. 

            Arguments:
                strands (List[DnaStrand]): The list of strands for the design.

            Strands may not have colors assigned to them. If they do then they have both
            an RGB and integer representation. The integer representation can be used
            as an ID to group staple strands by functionality.
        """
        for strand in strands:
            if (strand.is_scaffold):
                continue 
            base = strand.tour[0]
            for staple_color in self.staple_colors:
                if ((staple_color.vhelix_num == base.h) and (staple_color.vhelix_pos == base.p)): 
                    strand.color = staple_color.rgb 
                    strand.icolor = staple_color.color
            #__for color in self.staple_colors
        #__for strand in strands

    def set_sequence_from_name(self, dna_structure, modified_structure, seq_name):
        """ Set the sequence information for the staple and scaffold strands using a known
            origami vector sequence name.

            Cadnano seems to start assigning the sequence at virtual helix 0. If ordered_traverse=True
            then use the cadnano method.

            If the structure has not been modified with insertions and deletions then we will need to 
            use its insertions and deletions information (at the base level) to selectively set its 
            sequence.

            Arguments:
                modified_structure (bool): If true then the structure has been modified for insertions and deletions. 
                seq_name (string): The name of the sequence as defined in the dna_sequence_data dictionary. 
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.setLevel(logging.INFO)
        sequence = dna_sequence_data.get(seq_name,None)
        seq_index = 0
        sequence_length = len(sequence)
        self._logger.debug("-------------------- set_sequence_from_name p --------------------")
        self._logger.debug("sequence name: %s" % seq_name)
        self._logger.debug("sequence length: %d" % len(sequence))
        base_connectivity = dna_structure.base_connectivity
        self._logger.debug("Size of base_connectivity %d" % len(base_connectivity))
        strands = dna_structure.strands 
        ordered_traverse = True
        ordered_traverse = False

        for strand in strands:
            if (not strand.is_scaffold):
                continue 

            if (ordered_traverse): 
                self._logger.debug("---------- ordered traverse scaffold strand --------------------")
                tour = strand.tour
                base = tour[0]
                min_vh = base.h
                min_p = base.p
                start_index = 0

                for j in xrange(0,len(tour)):
                    base = tour[j]
                    if (base.h < min_vh):
                        min_vh = base.h
                        min_p = base.p
                        start_index = j
                    if ((base.h == min_vh) and (base.p < min_p)):
                        min_p = base.p
                        start_index = j
                self._logger.debug("start_index %d  min_vh %d  min_p %d ", start_index, min_vh,min_p)
                start_id = tour[start_index].id
                start_base = tour[start_index]
                base = start_base 
                even_vh = min_vh % 2

                for j in xrange(0,len(tour)):
                    letter = sequence[seq_index]
                    up_base = base.up
                    down_base = base.down
                    across_base = base.across

                    if (not modified_structure):
                        if (across_base != None):
                            if (base.num_deletions  != 0):
                                self._logger.debug("deleted base: id %d " % base.id)
                                letter = 'N'
                            elif (base.num_insertions != 0):
                                self._logger.debug("inserted base: id %d " % base.id)
                                strand.insert_seq.append(sequence[seq_index+1])
                                seq_index += 2
                            else:
                                seq_index += 1
                    else:
                        seq_index += 1

                    if (seq_index == sequence_length): 
                        seq_index = 0

                    base.seq = letter 
                    self._logger.debug("base id %d  vh %d  pos %d  up %d  down %d  across %d  seq %s", 
                        base.id, base.h, base.p, up, down, across, base.seq)

                    if (across_base != None):
                        across_base.seq = self._wspair(letter)

                    if (even_vh):
                        if (up != -1):
                            base = base_up
                        elif (across_base != None):
                            base = across_base 
                            even_vh = base.h % 2
                    else:
                        if (down_base != None):
                            base = down_base
                        elif (across_base != None):
                            base = across_base
                            even_vh = base.h % 2

                 #__for j in xrange(0,len(tour)):

            else: 
                self._logger.debug("---------- unordered traverse scaffold strand --------------------")
                for i in xrange(0,len(strand.tour)):
                    letter = sequence[seq_index]
                    base = strand.tour[i]
                    across_base = base.across

                    if (not modified_structure):
                        if (base.num_deletions != 0):
                            self._logger.debug("deleted base: id %d  vh %d  pos %d", base.id, base.h, base.p)
                            letter = 'N'
                        elif (base.num_insertions != 0):
                            self._logger.debug("inserted base: id %d  vh %d  pos %d", base.id, base.h, base.p)
                            strand.insert_seq.append(sequence[seq_index+1])
                            seq_index += 2
                        else:
                            seq_index += 1
                    else:
                        seq_index += 1

                    base.seq = letter 
                    self._logger.debug("base id %d  vh %d  pos %d  seq %s", base.id, base.h, base.p, base.seq)

                    if (seq_index == sequence_length): 
                        seq_index = 0

                    if (base.across != None):
                        base.across.seq = self._wspair(letter)
                #__for i in xrange(0,len(strand.tour))
        #__for strand in strands

        print_strands = True
        print_strands = False
        if (print_strands):
            self._logger.debug("---------- strands sequences ----------")
            self._logger.debug(">>> number of strands: %d " % len(strands))
            for strand in strands:
                tour = strand.tour 
                self._logger.debug(">>> strand: %d  scaf: %d length: %d" % (strand.id,strand.is_scaffold,len(tour)))
                self._logger.debug("    seq:")
                for base in tour:
                    self._logger.debug("    vhelix: %d  pos: %d  seq: %s" % (int(base.h), int(base.p), base.seq))
            #__for i in xrange(0,len(strands))

    #__def set_sequence_from_name

    def set_sequence(self, dna_structure, modified_structure, sequence):
        """ Set the sequence information for the staple and scaffold strands.

            The caDNAno csv file contains sequence information reflecting the insertions and deletions of a
            design. If a structure has had bases inserted and deleted then its sequence can be set directly 
            from the sequences in the csv file. If the structure has not been modified with insertions and
            deletions then we will need to use its insertions and deletions information (at the base level)
            to set selectively set its sequence from the sequences in the csv file.   

            Arguments:
               modified_structure (bool): If True then the structure has been modified with deletions and insertions. 
               sequence (DnaSequence): A list of DnaSequence objects representing the sequences for staple or scaffold strands.
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.setLevel(logging.INFO)

        strands = dna_structure.strands 
        staple_ends = dna_structure.staple_ends 

        for i in xrange(0,len(sequence)):
            seq = sequence[i]
            h0 = int(sequence[i].start[0])
            p0 = int(sequence[i].start[1])
            istrand = staple_ends[(h0,p0)]
            strand = strands[istrand-1]
            tour = strand.tour

            if (modified_structure):
                for j in xrange(0,len(strands[istrand-1].tour)):
                    base = tour[j]
                    base.seq = seq.letters[j]
                    if (base.across != None):
                        base.across.seq = self._wspair(seq.letters[j])
                #__for j

            else:
                seq_index = 0
                for j in xrange(0,len(tour)):
                    letter = seq.letters[seq_index]
                    base = tour[j]
                    if (base.num_deletions != 0):
                        letter = 'N'
                    elif (base.num_insertions != 0):
                        strand.insert_seq.append(seq.letters[seq_index+1])
                        seq_index += 2
                    else:
                        seq_index += 1
                    base.seq = letter
                    if (base.across != None):
                        base.across.seq = self._wspair(letter)
                #__for j

        #__for i
    #__def set_sequence

    def _wspair(self, x):
        """ Match a base with its complementary base. 
        """
        x = x.upper()

        if (x == 'A'):
            y = 'T'
        elif (x == 'G'):
            y = 'C'
        elif (x == 'C'):
            y = 'G'
        elif (x == 'T'):
            y = 'A'
        elif (x == 'N'):
            y = 'N'
        else:
            self._logger.error('Illegal base.')
        return y

    class StapleColor(object):
        """ This class stores data for a staple color.
        """
        def __init__(self,staple_color, vhelix_num):
            self.vhelix_num = vhelix_num 
            self.vhelix_pos = int(staple_color[0])
            self.color = staple_color[1]
            red = (self.color>>16) & 0xFF
            green = (self.color>>8) & 0xFF
            blue  = (self.color) & 0xFF
            self.rgb = [red/255.0, green/255.0, blue/255.0]
            #print("[DnaStapleColor] red %d  green %d  blue %d " % (red, green, blue))

        def match(self,pos):
            if (self.vhelix_pos == pos):
               return list(self.rgb)
            else:
               return []

#__class CadnanoTopology(object)

def main():
    """ Create a topology table from a caDNAno JSON design file."""
    json_file_name = sys.argv[1]
    cadnano_reader = CadnanoReader()
    #cadnano_reader.set_logging_level(logging.DEBUG)
    cadnano_model = cadnano_reader.read_json(json_file_name)
    convert_design = CadnanoConvertDesign()
    dna_structure = convert_design.create_structure(cadnano_model)

    if (len(sys.argv) == 3):
       csv_file_name = sys.argv[2]
       seq = cadnano_reader.read_csv(csv_file_name)
       structure.set_sequence(seq)

if __name__ == '__main__':
    main()

