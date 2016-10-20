#!/usr/bin/env python
""" 
This module is used to create a DNA structure from a caDNAno DNA origami design file. 

The virtual helices from the CadnanoDesign object are used to create a list of DnaStructureHelix objects. 
DnaBase objects are created from the scaffold and staple bases defined for each virtual helix and stored 
in the appropriate DnaStructureHelix object. A DnaStructure object object is created for the design. It 
stores the list of DnaStructureHelix objects and a list of all the bases defined for the design. 

This code is based on a direct translation of the set MATLAB scripts to convert a caDNAno design to a CanDo 
.cndo file from Mark Bathe's Laboratory for Computational Biology & Biophysics at MIT.
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
from math import sqrt,cos,acos,sin,asin,pi
from dna_sequence_data import dna_sequence_data

from design import CadnanoDesign,CadnanoVirtualHelix,CadnanoBase
from reader import CadnanoReader 
from common import CadnanoLatticeType

try:
    import os.path
    base_path = os.path.abspath( os.path.dirname(__file__) + '/../../../' )
    sys.path.append( base_path )

    from nanodesign_transition.base import DnaBase 
    from nanodesign_transition.strand import DnaStrand
    from nanodesign_transition.dna_structure import DnaStructure,DnaStructureHelix
    from nanodesign_transition.lattice import Lattice,SquareLattice,HoneycombLattice
    from nanodesign_transition.parameters import DnaPolarity,DnaParameters

    sys.path = sys.path[:-1]
except ImportError as i:
    print "Could not get nanodesign_transition module"
    raise i

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
        self._timer = _Timer()
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
        console_handler = logging.StreamHandler()
        #formatter = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s - %(message)s')
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
        self._logger.info("Distance between adjacent helices %g " % self.dna_parameters.helix_distance)
        self._logger.info("Helix radius %g " % self.dna_parameters.helix_radius) 
        self._timer.start()

        # Create a list of DnaStructureHelix objects for the design. 
        helices = self._create_structure_topology_and_geometry(design)
        self._logger.info("Number of bases in design %d " % len(self.base_map))
        self._logger.info("Time to create structure topology table %s " % self._timer.finish())

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
        print_base_connectivity = False
        if print_base_connectivity:
            self._logger.info("---------- base_connectivity  ---------- ")
            self._logger.info("size of base_connectivity: %d " % len(base_connectivity))
            for i in xrange(0,num_bases):
                base = base_connectivity[i]
                up = base.up.id if base.up else -1
                down = base.down.id if base.down else -1
                across = base.across.id if base.across else -1
                self._logger.info("%4d: base id %4d  up %4s  down %4s  across %4s  h %4d  p %4d  scaf %d" %
                    ( i, base.id, up, down, across, base.h, base.p, base.is_scaf))
        #__if print_base_connectivity_p

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
        base_connectivity,strands = self._build_strands(base_connectivity)
        self._set_strands_colors(strands)
        self.dna_structure.strands = strands
        self._logger.info("Number of strands %d " % len(strands)) 
        print_strands = False
        if print_strands:
            for strand in strands:
                self._logger.info("Strand %4d  bases %4d  scaf %6s  start helix %4d  pos %4d" % (strand.id, 
                    len(strand.tour), strand.is_scaffold, strand.tour[0].h,  strand.tour[0].p))
        #__if print_strands

        # Calculate staple ends.
        staple_ends = self._calculate_staple_ends(base_connectivity)
        self.dna_structure.staple_ends = staple_ends

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
            for base in helix.staple_base_list:
                if base.num_insertions != 0:
                    num_insert_bases += base.num_insertions 
                    base_inserts[base.id] = base
            #__for base in helix.staple_base_list
            for base in helix.scaffold_base_list:
                if base.num_insertions != 0:
                    base_inserts[base.id] = base
                    self._logger.debug("Insert base id %d  h %d  pos %d" % (base.id, base.h, base.p))
            #__for base in helix.scaffold_base_list
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
        rot_mat = _vrrotvec2mat(y_up_vec, _deg2rad(ang_bp))
        next_frame = np.dot(rot_mat,curr_frame)
        [insert_coords, insert_frames] = _bp_interp(curr_coords, curr_frame, next_coords, next_frame, num_inserts)
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
        rot_mat = _vrrotvec2mat(y_up_vec, _deg2rad(ang_bp))
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
            if (False):
                s = ""
                for base in scaffold_bases:
                    s += str(base.p) + " "
                self._logger.debug("Scaffold bases positions %s " % s) 

	    self._logger.debug("Number of staple bases %d " % len(staple_bases)) 
            if (False):
                s = ""
                for base in staple_bases:
                    s += str(base.p) + " "
                self._logger.debug("Staple bases positions %s " % s) 
   
            # Generate the helix axis coordinates and frames, and DNA helix nucleotide coordinates.
            axis_coords, axis_frames, scaffold_coords, staple_coords = \
                self._generate_coordinates(lattice_type, row, col, num, scaffold_bases, staple_bases)

            # Create a dna structure object that stores the helix information. 
            structure_helix = DnaStructureHelix(i, num, scaffold_polarity, axis_coords, axis_frames, scaffold_bases, 
                                                staple_bases)
            structure_helix.lattice_num = num
            structure_helix.lattice_row = row
            structure_helix.lattice_col = col
            structure_helix.lattice_max_vhelix_size = max_vhelix_size 
            structure_helix.staple_base_list = staple_bases
            structure_helix.scaffold_base_list = scaffold_bases
            structure_helices.append(structure_helix)
        #__for vhelix in vhelices
        return structure_helices

    def _create_nt_map_table(self, base_connectivity, dnode, id_nt_0):
        """ Create the map table (id_nt). """
        self._logger.info("================================= _create_nt_map_table =======================")
        n_bp = dnode.shape[0]
        id_nt = np.zeros((n_bp, 2), dtype=int)
        self._logger.info("n_bp %d" % n_bp)

        num_pairs = 0
        for i,base in enumerate(base_connectivity):
           if not (base.across and base.is_scaf): 
               continue 
           id_nt[num_pairs,0] = base.id 
           id_nt[num_pairs,1] = base.across.id 
           num_pairs += 1
        self._logger.info("num_pairs %d" % num_pairs)
        return id_nt 

        for i in xrange(0,n_bp):
            # Scaffold nucleotide
            base_data = id_nt_0[i,0:3]
            id_nt[i,0] = self._find_base_location(base_data, base_maps)
            # Staple nucleotide
            base_data = id_nt_0[i,3:6]
            id_nt[i,1] = self._find_base_location(base_data, base_maps)
        return id_nt 

    def _find_base_location(self, base_data, base_maps):
        """ Find the location in the base topology array of the given base using the helix number it is in 
            and its position within that helix. 

            This is used to convert a helix-based base indexing using a helix number and position to a global ID.
            The query is a tuple (helix num, position). For a base with no up, down or across base the query
            will be (-1,-1) and the location returned -1.

            Arguments:
                base_data (list(float)[3]): A list of three values representing the helix number, position 
                    within that helix and a flag if it is from a scaffold(0) or a staple(1) strand.
                base_maps (list(dict)[2]): A list of two base maps: 0:scaffold, 1:staple. 

            Returns the index into the base topology array of the given base.
        """
        query = tuple(int(base_data[j]) for j in xrange(0,2))
        strand_type = int(base_data[2])
        base_map = base_maps[strand_type]
        if query not in base_map:
            loc = -1
        else:
            loc = base_map[query]
        return loc

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

    def _generate_coordinates(self, lattice_type, row, col, helix_num, scaffold_bases, staple_bases):
        """ Generate the axis coordinates, axis reference frames, and nucleotide coordinate for the
            virtual helix. 

            Arguments:
                lattice_type (CadnanoLatticeType): The lattice type for this design.
                row (int): The caDNAno row number. 
                col (int): The caDNAno column number. 
                helix_num (int): The caDNAno virtual helix number. 
                scaffold_bases (List[DnaBase]): The list of scaffold bases for this helix. 
                staple_bases (List[DnaBase]): The list of staple bases for this helix. 

            The helix coordinates and reference frames are generated for all of the virtural helix
            positions that contain a base. The coordinates and refereance frames are also set for 
            scaffold and staple bases. 
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.setLevel(logging.INFO)
        self._logger.debug("-------------------- generate_coordinates --------------------")
        r_strand = self.dna_parameters.helix_distance / 2.0 # half the distance between the axes of adjacent DNA helices.
        r_helix = self.dna_parameters.helix_radius          # radius of DNA helices (nm)
        dist_bp = self.dna_parameters.base_pair_rise        # rise between two neighboring base-pairs (nm)
        ang_bp = self.dna_parameters.base_pair_twist_angle  # twisting angle between two neighboring base-pairs (degrees)
        ang_minor = self.dna_parameters.minor_groove_angle  # angle of the minor groove (degrees)

        # Positions of the scaffold nucleotide and staple nucleotide
        # in the local reference frame.
        scaf_local = r_helix * np.array([cos(_deg2rad(180-ang_minor/2)), sin(_deg2rad(180-ang_minor/2)), 0.0]).transpose()
        stap_local = r_helix * np.array([cos(_deg2rad(180+ang_minor/2)), sin(_deg2rad(180+ang_minor/2)), 0.0]).transpose()

        # Set the helix start coordinates, based on helix (row,col), and frame orientation, 
        # based on helix number (polarity).
        if (lattice_type == CadnanoLatticeType.honeycomb):
            xpos =  sqrt(3.0) * col * r_strand
            zpos = -3.0 * row * r_strand;
            if ( ((row % 2 == 0) and (col % 2 == 0)) or ((row % 2 != 0) and (col % 2 != 0))):
                zpos = zpos + r_strand
            if ( helix_num % 2 == 0 ):
                init_ang = -30 + ang_bp / 2
            else:
                init_ang = 150 + ang_bp / 2
        elif (lattice_type == CadnanoLatticeType.square): 
            xpos =  2.0 * col * r_strand
            zpos = -2.0 * row * r_strand
            if helix_num % 2 == 0:
                init_ang = 180 + ang_bp / 2
            else:
                init_ang = 0 + ang_bp / 2
        #__if (lattice_type == CadnanoLatticeType.honeycomb)

        init_coord_ang_strand = np.array([xpos, 0.0, zpos, init_ang]);

        if (helix_num % 2 == 0):
            e3 = np.array([0, 1, 0],dtype=float)
        else:
            e3 = np.array([0, -1, 0],dtype=float)

        # Create a list of sorted base positions.
        base_positions = set()
        for base in scaffold_bases: 
            base_positions.add(base.p)
        for base in staple_bases: 
            base_positions.add(base.p)
        #self._logger.debug("Base positions %s" % str(sorted_base_positions))

        # Compute helix axis coordinates and frames.
        num_base_positions = len(base_positions)
        axis_coords = np.zeros((num_base_positions,3), dtype=float)
        axis_frames = np.zeros((3,3,num_base_positions), dtype=float);
        pos_map = {}
        self._logger.debug("Axis frames: ")
        for i,p in enumerate(sorted(base_positions)): 
            pos_map[p] = i
            ref_coord = init_coord_ang_strand + np.array([0, dist_bp*p, 0, ang_bp*p]);
            # Base coordinate. 
            axis_coords[i,:] = ref_coord[0:3];
            # Base orientation.
            e2 = np.array([cos(-_deg2rad(ref_coord[3])), 0, sin(-_deg2rad(ref_coord[3]))])
            e1 = np.cross(e2, e3);
            axis_frames[:,:,i] = np.array([e1, e2, e3]).transpose();
            self._logger.debug("pos %d  frame 1 %g %g %g " % (p, axis_frames[0,0,i], axis_frames[1,0,i], axis_frames[2,0,i]))
            self._logger.debug("        frame 2 %g %g %g " % (   axis_frames[0,1,i], axis_frames[1,1,i], axis_frames[2,1,i]))
            self._logger.debug("        frame 3 %g %g %g " % (   axis_frames[0,2,i], axis_frames[1,2,i], axis_frames[2,2,i]))
        #__for base in staple_bases

        # Compute scaffold nucleotide positions and set base coordinates and frame.
        num_scaffold_bases = len(scaffold_bases)
        scaffold_coords = np.zeros((num_scaffold_bases,3), dtype=float)
        for i,base in enumerate(scaffold_bases): 
            j = pos_map[base.p]
            base.coordinates = axis_coords[j]
            base.ref_frame = axis_frames[:,:,j]
            scaffold_coords[i,:] = axis_coords[j,:] + np.dot(axis_frames[:,:,j], scaf_local)
            base.nt_coords = scaffold_coords[i]
        #__for base in scaffold_bases 

        # Compute staple nucleotide positions and set base coordinates and frame.
        num_staple_bases = len(staple_bases)
        staple_coords = np.zeros((num_staple_bases,3), dtype=float)
        for i,base in enumerate(staple_bases): 
            j = pos_map[base.p]
            base.coordinates = axis_coords[j]
            base.ref_frame = axis_frames[:,:,j]
            staple_coords[i,:] = axis_coords[j,:] + np.dot(axis_frames[:,:,j], stap_local)
            base.nt_coords = staple_coords[i]
        #__for base in stape_bases

        return axis_coords, axis_frames, scaffold_coords, staple_coords 

    def _calculate_lattice_directions(self):
        """ Calculate the unit vectors pointing to neighboring cells.
        """
        if (lattice_type == CadnanoLatticeType.honeycomb):
            dx =  sqrt(3.0) 
            dz = -3.0 

        elif (lattice_type == CadnanoLatticeType.square):
            dx =  2.0 
            dz = -2.0 

    def _calculate_staple_ends(self, base_connectivity):
        base_connectivity,strands = self._build_strands(base_connectivity)
        num_strands = len(strands)
        staple_ends = np.empty(shape=(0,5),dtype=float)

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
                staple_ends = np.concatenate((staple_ends, [[i+1, h0, p0, h1, p1]]), axis=0)
        #__for i in xrange(0,num_strands)__
        return staple_ends

    def _set_possible_crossovers(self,design):
        """ Set the possible cross-overs for scaffold and staple strands.
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.setLevel(logging.INFO)
        self._logger.debug("-------------------- set_possible_crossovers --------------------")
        scaffold_crossovers = []
        staple_crossovers = [] 
        #structure_helices = self.dna_structure.structure_helices
        structure_helices_coord_map = self.dna_structure.structure_helices_coord_map
        for vhelix in design.helices:
            num = vhelix.num 
            col = vhelix.col
            row = vhelix.row
            self._logger.debug(">>> vhelix: num: %d  row: %d  col: %d " % (num, row, col))
            staple_crossovers = vhelix.possible_staple_crossovers
            scaffold_crossovers = vhelix.possible_scaffold_crossovers
            self._logger.debug("            num staple cross-overs: %d " % len(staple_crossovers)) 
            shelix = structure_helices_coord_map[(row,col)]
            for cross_vh,index in staple_crossovers:
                self._logger.debug("            staple cross-over: %d,%d " % (cross_vh.num,index))
                cross_sh = structure_helices_coord_map[(cross_vh.row,cross_vh.col)]
                shelix.possible_staple_crossovers.append((cross_sh,index))
            self._logger.debug("            num scaffold cross-overs: %d " % len(scaffold_crossovers)) 
            for cross_vh,index in scaffold_crossovers:
                self._logger.debug("            scaffold cross-over: %d,%d " % (cross_vh.num,index))
                cross_sh = structure_helices_coord_map[(cross_vh.row,cross_vh.col)]
                shelix.possible_scaffold_crossovers.append((cross_sh,index))
        #__for vhelix in design.helices

    def _build_strands(self, base_connectivity):
        #self._logger.setLevel(logging.DEBUG)
        self._logger.setLevel(logging.INFO)
        self._logger.debug("==================== build strands ====================")
        num_bases = len(base_connectivity)
        self._logger.debug("Number of bases %d" % num_bases) 
        strands = []
        n_strand = 0
        is_visited = [False]*num_bases

        while (True):
            try:
                base_index = is_visited.index(False)
                curr_base = base_connectivity[base_index]
                self._logger.debug("------------- find first base --------------")
                self._logger.debug("Base_index %d" % base_index)
                self._logger.debug("Curr_base id %d  h %d  p %d  up %s" % (curr_base.id, curr_base.h, curr_base.p, curr_base.up))
            except:
                break

            init_base = curr_base
            #init_baseID = curr_base.id

            # Find the first base in the current strand.
            self._logger.debug("------------- walk bases --------------")
            while curr_base.up and (curr_base.up.id != init_base.id):
                curr_base = curr_base.up
                self._logger.debug(" curr_base id %d  h %d  p %d" % (curr_base.id, curr_base.h,  curr_base.p))
                if (is_visited[curr_base.id]):
                    print('[build_strand] **** ERROR: Reached a visited base.');
                    return None,None
            #_while 
    
            self._logger.debug("---------- add strand %d ----------" % n_strand)
            self._logger.debug("first strand base: curr_base.id %d" % curr_base.id)
            strand = DnaStrand(n_strand,self.dna_structure)
            strands.append(strand)
    
            if not curr_base.up:                     # currBase is at the 5'-end of the strand
                strand.is_circular = False
            elif curr_base.up.id == init_base.id:    # currBase goes back to the starting point
                strand.is_circular = True
                curr_base = init_base
                self._logger.debug("strand is circular. " )
                self._logger.debug("curr_base.id %d" % curr_base.id)
            else:
                print('[build_strand] **** ERROR: Exception.')
                return None,None
    
            # Walk through the current strand.
            n_residue = 1
            strand.tour.append(curr_base)
            curr_base.strand = n_strand
            curr_base.residue = n_residue
            is_visited[curr_base.id] = True
    
            while ( (not strand.is_circular and curr_base.down) or 
                    (strand.is_circular and (curr_base.down.id != init_base.id)) ):
                curr_base = curr_base.down
                #if debug: print "[build_strand] curr_base id %d  h %d  p %d" % (curr_base.id, curr_base.h, curr_base.p )
    
                if (is_visited[curr_base.id]):
                    sys.stderr.write('[build_strand] **** ERROR: Reached a visited base.\n')
                    return None,None
    
                if (n_residue == 1):
                    strand.is_scaffold = curr_base.is_scaf
    
                n_residue = n_residue + 1
                strand.tour.append(curr_base)
                curr_base.strand = n_strand
                curr_base.residue = n_residue
                is_visited[curr_base.id] = True
            #__while((not strand.is_circular 

            # Modify the strand if it is circular and the first base crosses over to another helix.
            # The first base is moved to the end of the strand. This will prevent later issues, like
            # domains that don't follow the strand tour because the strand first domain would be
            # merged with the last domain during domain calculation.
            if strand.is_circular and len(strand.tour) > 2:
                first_base = strand.tour[0]
                second_base = strand.tour[1]
                if first_base.h != second_base.h:
                    del(strand.tour[0])
                    strand.tour.append(first_base)
            #__if strand.is_circular and len(strand.tour) > 2

            n_strand += 1
        #__while (True):
        return base_connectivity,strands

    def _set_strands_colors(self, strands):
        """ Set the color for staple strands. 

            Arguments:
                strands (List[DnaStrand]): The list of strands for the design.
        """
        for strand in strands:
            if (strand.is_scaffold):
                continue 
            base = strand.tour[0]
            for color in self.staple_colors:
                if ((color.vhelix_num == base.h) and (color.vhelix_pos == base.p)): 
                    strand.color = color.rgb 
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
        start = np.array([0,0], dtype=float)

        for i in xrange(0,len(sequence)):
            seq = sequence[i]
            start[0] = sequence[i].start[0]
            start[1] = sequence[i].start[1]
            row = _find_row(start, staple_ends[:,1:3])[0]
            istrand = int(staple_ends[row,0])
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

        print_strands = True
        print_strands = False
        if (print_strands):
            self._logger.debug("---------- strands sequences ----------")
            self._logger.debug(">>> number of strands: %d " % len(strands))
            for strand in strands:
                tour = strand.tour
                self._logger.debug(">>> strand: %d scaf: %d len: %d" % (strand.id,strand.is_scaffold,len(tour)))
                self._logger.debug("    seq:")
                for base in tour:
                    self._logger.debug("    vhelix: %d  pos: %d  seq: %s" % (int(base.h), int(base.p), base.seq))
            #__for i in xrange(0,len(strands))
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




class _Timer(object):
    """ The Timer class is used to calculare elapsed time between calls to the start and finish methods."""
    def __init__(self):
        self.start_time = 0.0
        self.end_time = 0
        self.secs = 0
        self.msecs = 0

    def start(self):
        self.start_time = time.time()

    def finish(self):
        self.end_time = time.time()
        self.secs = self.end_time - self.start_time
        self.msecs = self.secs * 1000  # millisecs
        return self.secs

def _deg2rad(deg):
    """Convert degrees into radians. """
    rad = (pi/180)* deg;
    return rad

def _find_row(neigh, curr_bases):
    tmp = curr_bases - neigh
    s = np.sum(np.abs(tmp),1)
    ind = np.where(s==0)[0]
    return ind



def _vrrotmat2vec(R):
    """ Extract the equivalent rotation about an axis from a rotation matrix. """
    m00 = R[0,0]
    m01 = R[0,1]
    m02 = R[0,2]
    m10 = R[1,0]
    m11 = R[1,1]
    m12 = R[1,2]
    m20 = R[2,0]
    m21 = R[2,1]
    m22 = R[2,2]
    angle = acos(( m00 + m11 + m22 - 1)/2.0)
    x = (m21 - m12) / sqrt( pow(m21-m12,2) + pow(m02-m20,2) + pow(m10-m01,2) )
    y = (m02 - m20) / sqrt( pow(m21-m12,2) + pow(m02-m20,2) + pow(m10-m01,2) )
    z = (m10 - m01) / sqrt( pow(m21-m12,2) + pow(m02-m20,2) + pow(m10-m01,2) )
    return np.array([x,y,z],dtype=float),angle


def _vrrotvec2mat(axis, theta):
    """ Create a rotation matrix to rotate theta degrees about the axis defined by vec. """
    s = np.sin(theta)
    c = np.cos(theta)
    t = 1 - c
    #print("[vrrotvec2mat] theta=%s" % str(theta))

    # normalize the vector
    x,y,z = axis / np.linalg.norm(axis)
    #print("[vrrotvec2mat] x=%s" % str(x))
    #print("[vrrotvec2mat] y=%s" % str(y))
    #print("[vrrotvec2mat] z=%s" % str(z))

    return np.array([ [t*x*x + c,   t*x*y - s*z,  t*x*z + s*y],
                      [t*x*y + s*z, t*y*y + c,    t*y*z - s*x],
                      [t*x*z - s*y, t*y*z + s*x,  t*z*z + c  ]])


def _bp_interp(dnode_1, triad_1, dnode_2, triad_2, n):
    """ Interpolate the position (dnode) and orientation (triad) between two base pairs.

        Solve for rotation matrix R defined as the solution of: 
            R * triad_1 = triad_2 (multiply both sides by transpose(triad_1).
    """
    dnode_interp = np.zeros((n,3),dtype=float)
    triad_interp = np.zeros((3,3,n),dtype=float)

    # solve for rotation matrix R
    R = np.dot(triad_2,triad_1.T)
    #print("[bp_interp] R=%s " % str(R))
    #print("[bp_interp] R.shape=%s " % str(R.shape))

    a,theta = _vrrotmat2vec(R)
    #print("[bp_interp] a=%s " % str(a))
    #print("[bp_interp] theta=%s " % str(theta))

    # Calculate for dNode_interp and triad_interp
    for i in xrange(0,n):
        dnode_interp[i,:] = (dnode_1*(n+1-(i+1)) + dnode_2*(i+1)) / (n+1)
        angle = theta*(i+1)/(n+1)
        #print("[bp_interp] angle=%s " % str(angle))
        rot_mat = _vrrotvec2mat(a, angle)
        #print("[bp_interp] rot_mat=%s\n" % str(rot_mat))
        #print("[bp_interp] rot_mat.shape=%s\n" % str(rot_mat.shape))
        triad_interp[:,:,i] = np.dot(rot_mat,triad_1)

    return dnode_interp, triad_interp


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

