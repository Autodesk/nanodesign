# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This module defines the classes used to define a structure helix of a DNA structure.

A structure helix is a region in a DNA structure that forms a cylindrical structural element. It can be composed of 
one or two DNA strands defined by the set of scaffold and staple bases assigned to the helix. 

A helix can be considered to be a container for storing scaffold and staple bases. The position of a base in a helix 
is given by an integer between 0 and N-1, where N is the maximum number of bases staple or scaffold, assigned to it. 

    Positions     0       1       2       3       4       5             N-1
              -----------------------------------------------------------------
    scaffold  | base1 | base2 | base3 | base4 | base5 | base6 |  ...  | baseN |
              -----------------------------------------------------------------
    staple    | base1 | base2 | base3 | base4 | base5 | base6 |       | baseN |
              -----------------------------------------------------------------
"""
from collections import OrderedDict
import inspect
import itertools
import json
import logging
import numpy as np
import sys

# Within package imports
from .parameters import DnaParameters,DnaPolarity
from ..converters.cadnano.common import CadnanoLatticeType
from .lattice import Lattice
from .base import DnaBase

class DnaStructureHelix(object):
    """ This class stores information for a DNA structure helix. 

        Attributes:
            load_order (int): The order (0-based) the helix was processed. This is needed when writing a caDNAno file. 
            dna_parameters (DnaParameters): The DNA parameters to use when creating the 3D geometry for the design.
            end_coordinates (NumPy 2x3 ndarray[float]): The coordinates at the ends of the helix.
            end_frames (NumPy 3x3x2 ndarray[float]): The coordinate frames at the ends of the helix.
            helix_axis_frames (NumPy 3x3xN ndarray[float]): The coordinate frames of base nodes along the helix axis, 
                where N is the number of bases of the caDNAno virtual helix.
            helix_axis_coords (NumPy Nx3 ndarray[float]): The coordinates of base nodes along the helix axis, where N is the 
                number of bases of the caDNAno virtual helix. The elements of the array are sorted according to their position
                in the helix. 
            helix_connectivity (List[DnaHelixConnection]): The list helices that this helix is connected to.
            id (int): Helix ID. This should be the same IDs used for a base helix number (i.e. the DnaBase.h ID).
            lattice_col (int): The caDNAno lattice column number.
            lattice_num (int): The caDNAno helix number.
            lattice_row (int): The caDNAno lattice row number. 
            lattice_type (CadnanoLatticeType): The lattice type the geometry of this structure is derived from.
            possible_scaffold_crossovers (List[(DnaStructureHelix,int,NumPy 3x1 array[float])]): The list of possible 
                scaffold crossovers for this helix. Each entry in the list is a 3-tuple giving the helix it crosses over 
                to (object pointer), the helix position where the crossover occurs and the helix axis coodinates.
            possible_staple_crossovers (List[(DnaStructureHelix,int,NumPy 3x1 array[float])]): The list of possible staple 
                crossovers for this helix. Each entry in the list is a 3-tuple giving the helix it crosses over to 
                (object pointer), the helix position where the crossover occurs the helix axis coodinates.
            scaffold_bases (List[DnaBase]): The list of helix scaffold bases. 
            scaffold_polarity (DnaPolarity): The polarity of the scaffold. 
            scaffold_pos (Dict[int,DnaBase]: The dict mapping helix positions to scaffold bases.
            staple_bases (List[DnaBase]): The list of helix staple bases. 
            staple_pos (Dict[int,DnaBase]: The dict mapping helix positions to staple bases.

        The geometry of a helix is determined by the coordinates of base nodes along its axis. The coordinates and reference
        frames of the bases defined for the helix are references into the helix_axis_coords and helix_axis_frames arrays.
    """ 

    # TODO (DaveP) We need to remove the references to caDNAno lattice-based information.

    def __init__(self, load_order, id, scaffold_polarity, helix_axis_coords, helix_axis_frames, scaffold_coords, staple_coords,
                 scaffold_bases, staple_bases):
        """ Initialize a DnaStructureHelix object.

            Arguments:
                load_order (int): The order the helix was processed.
                id (int): The helix ID.
                scaffold_polarity (DnaPolarity): The helix polarity. 
                helix_axis_coords (NumPy Nx3 ndarray[float]): The coordinates of base nodes along the helix axis. 
                helix_axis_frames (NumPy 3x3xN ndarray[float]): The coordinate frames base nodes alonge the helix axis.
                scaffold_bases (List[DnaBase]): The list of scaffold bases defined for the helix.
                scaffold_coords (NumPy Nx3 ndarray[float]): The scaffold DNA helix nucleotide coordinates. 
                staple_bases (List[DnaBase]): The list of staple bases defined for the helix.
                staple_coords (NumPy Nx3 ndarray[float]): The staple DNA helix nucleotide coordinates. 
        """
        self.id = id
        self.load_order = load_order
        self.staple_bases = staple_bases
        self.staple_coords = staple_coords
        self.scaffold_bases = scaffold_bases
        self.scaffold_coords = scaffold_coords
        self.helix_axis_frames = helix_axis_frames
        self.helix_axis_coords = helix_axis_coords
        self.scaffold_polarity = scaffold_polarity 
        self.lattice_row = -1
        self.lattice_col = -1
        self.lattice_num = -1
        self.lattice_max_vhelix_size = 0
        self.helix_connectivity = []
        self.possible_staple_crossovers = []
        self.possible_scaffold_crossovers = []
        self._logger = logging.getLogger(__name__ + ":" + str(self.id))

        # Set helix ends coordinates.
        self.end_coordinates = np.zeros((2,3), dtype=float)
        self.end_frames = np.zeros((3,3,2), dtype=float)
        self.staple_pos = {}
        self.scaffold_pos = {}
        self.set_end_coords()
        self.build_base_pos_maps()

    def set_end_coords(self):
        """ Set the helix end coordinates and reference frames. """

        if len(self.scaffold_bases) != 0:
            point1 = self.scaffold_bases[0].coordinates
            point2 = self.scaffold_bases[-1].coordinates
            frame1 = self.scaffold_bases[0].ref_frame
            frame2 = self.scaffold_bases[-1].ref_frame
        else:
            point1 = self.staple_bases[0].coordinates
            point2 = self.staple_bases[-1].coordinates
            frame1 = self.staple_bases[0].ref_frame
            frame2 = self.staple_bases[-1].ref_frame

        for i in xrange(0,3):
            self.end_coordinates[0,i] = point1[i]
            self.end_coordinates[1,i] = point2[i]

        self.end_frames[:,:,0] = frame1
        self.end_frames[:,:,1] = frame2
    #__def set_end_coords

    def set_coordinates(self, helix_axis_coords, helix_axis_frames, scaffold_coords, staple_coords):
        """ Set the helix axis coordinates and frames, and DNA helix nucleotide coordinates. 

            Arguments:
                helix_axis_coords (NumPy Nx3 ndarray[float]): The coordinates of base nodes along the helix axis. 
                helix_axis_frames (NumPy 3x3xN ndarray[float]): The coordinate frames base nodes alonge the helix axis.
                scaffold_coords (NumPy Nx3 ndarray[float]): The scaffold DNA helix nucleotide coordinates. 
                staple_coords (NumPy Nx3 ndarray[float]): The staple DNA helix nucleotide coordinates. 
        """
        self.helix_axis_coords = helix_axis_coords
        self.helix_axis_frames = helix_axis_frames
        self.staple_coords = staple_coords
        self.scaffold_coords = scaffold_coords
        self.set_end_coords()
    #__def set_coordinates

    def get_start_pos(self):
        """ Get the starting helix position of the scaffold or staple strands. 
        """
        num_bases = len(self.staple_bases)
        if num_bases == 0:
            return None 
        staple_start_pos = self.staple_bases[0].p
        scaffold_start_pos = self.scaffold_bases[0].p
        start_pos = min(staple_start_pos, scaffold_start_pos)
        return start_pos 
    #__def get_start_pos

    def has_base_pos(self, pos):
        """ Check if the helix has bases at the given helix position.

            Arguments:
                pos (int): The position in the helix to check. 

            Returns two booleans, each set to True if there is a staple base and a scaffold 
            base at the given position.
        """
        return pos in self.staple_pos, pos in self.scaffold_pos
    #__def has_base_pos

    def build_base_pos_maps(self):
        """ Create maps for base positions. 

            This creates dicts to mapping helix positions to bases:
                self.staple_pos
                self.scaffold_pos

            The max/min positions for staple and scaffold bases are also calculated.
                self.min_staple_pos 
                self.max_staple_pos 
                self.min_scaffold_pos 
                self.max_scaffold_pos 
        """
        self.staple_pos = {}
        self.min_staple_pos = None
        self.max_staple_pos = None
        for base in self.staple_bases:
            self.staple_pos[base.p] = base
            if (self.min_staple_pos == None) or (base.p < self.min_staple_pos): 
                self.min_staple_pos = base.p
            if (self.max_staple_pos == None) or (base.p > self.max_staple_pos): 
                self.max_staple_pos = base.p
        #__for base in self.staple_bases
        self.scaffold_pos = {}
        self.min_scaffold_pos = None
        self.max_scaffold_pos = None
        for base in self.scaffold_bases:
            self.scaffold_pos[base.p] = base
            if (self.min_scaffold_pos == None) or (base.p < self.min_scaffold_pos): 
                self.min_scaffold_pos = base.p
            if (self.max_scaffold_pos == None) or (base.p > self.max_scaffold_pos): 
                self.max_scaffold_pos = base.p
        #__for base in self.scaffold_bases
    #__def build_base_pos_maps

    def add_maximal_staple_bases(self):
        """ Add bases for the maximal staple set. 

            New bases are created and added to the helix at locations not occupied by bases. 
            The helix locations to add bases are detemined by the max/min positions of the
            helix scaffold bases.
        """
        self._logger.debug("=================== add maximal staple bases %d ===================" % self.id)
        self._logger.debug("Scaffold polarity %s" % self.scaffold_polarity)
  
        # Determine the position to start adding bases.
        start_pos = self.min_scaffold_pos
        self._logger.debug("Start base pos %d " % start_pos) 

        # Determine the position to end adding bases
        end_pos = self.max_scaffold_pos
        self._logger.debug("End base pos %d " % end_pos) 

        # Add new bases.
        id = 0
        new_base_pos = set()
        for pos in xrange(start_pos,end_pos+1):
            has_staple,has_scaffold = self.has_base_pos(pos)
            if (not has_staple) and has_scaffold:
                base = DnaBase(id)
                base.p = pos
                base.h = self.id
                base.is_scaf = False
                base.across = self.scaffold_pos[pos]
                base.coordinates = np.array([0.0, 0.0, 0.0], dtype=float)
                self.staple_pos[base.p] = base
                new_base_pos.add(pos)
                id += 1
        #__for pos in xrange(start_pos,self.max_staple_pos)

        # Recreate list of staple bases.
        self.staple_bases = []
        for pos in sorted(self.staple_pos):
            base = self.staple_pos[pos]
            self.staple_bases.append(base)

        if self.scaffold_polarity == DnaPolarity.FIVE_PRIME:
            inc = 1
        else:
            inc = -1

        # Set base connectivity for new bases.
        for base in self.staple_bases:
            if base.p in new_base_pos: 
                up_pos = base.p + inc
                up_base = self.staple_pos.get(up_pos,None)
                if up_pos in new_base_pos:
                    base.up = up_base
                down_pos = base.p - inc
                down_base = self.staple_pos.get(down_pos,None)
                if down_pos in new_base_pos:
                    base.down = down_base
        #__for base in self.staple_bases

    #__def add_maximal_staple_bases

    def add_maximal_staple_crossovers(self):
        """ Add crossover connections for bases for the maximal staple set. """
        self._logger.debug("=================== add maximal staple crossovers %d ===================" % self.id)
        self._logger.debug("Scaffold polarity %s" % self.scaffold_polarity)

        # Create a dict mapping crossover positions to crossovers.
        crossover_pos = {}
        for crossover in self.possible_staple_crossovers:
            crossover_pos[crossover[1]] = crossover

        # Add the up/down pointers for bases at crossovers. 
        five_prime = (self.scaffold_polarity == DnaPolarity.FIVE_PRIME)
        for crossover in self.possible_staple_crossovers:
            to_helix = crossover[0]
            pos = crossover[1]
            has_staple,has_scaffold = self.has_base_pos(pos)
            if has_staple and has_scaffold and (pos in to_helix.staple_pos):
                base = self.staple_pos[pos]
                to_base = to_helix.staple_pos[pos]
                if five_prime:
                    if (base.p+1) in crossover_pos:
                        base.up = to_base 
                        to_base.down = base 
                    else:
                        base.down = to_base 
                        to_base.up = base 
                else:
                    if (base.p+1) in crossover_pos:
                        base.down = to_base 
                        to_base.up = base 
                    else:
                        base.up = to_base 
                        to_base.down = base 
                self._logger.debug("Crossover base at pos %d to helix %d" % (base.p, to_helix.id))
        #__for crossover in possible_staple_crossovers

    #__def add_maximal_staple_crossovers

    def apply_xform(self, xform):
        """ Apply a transformation to the helix coordiates and reference frames.

            Arguments:
                xform (Xform): The transformation to apply to the helx geometry.
        """
        self._logger.debug("=================== apply xform to helix %d ===================" % self.id)
        R = xform.rotation_matrix
        #xform.print_transformation()
        translation = xform.translation 
        center = xform.center 
        coord = np.array([0.0,0.0,0.0], dtype=float) 
        frame = np.zeros((3,3), dtype=float)
        xcoord = np.array([0.0,0.0,0.0], dtype=float) 
        xframe = np.zeros((3,3), dtype=float)
        num_coords = len(self.helix_axis_coords)
        self._logger.debug("Number of coordinates %d" % num_coords) 
        self._logger.debug("Xform center (%g %g %g)" % (center[0], center[1], center[2])) 
        self._logger.debug("Xform translation (%g %g %g)" % (translation[0], translation[1], translation[2])) 
        for i in xrange(0, num_coords):
            coord[:] = self.helix_axis_coords[i,:] - center
            frame = self.helix_axis_frames[:,:,i]
            xcoord[:] = np.dot(R, coord) + center + translation
            xframe = np.dot(R, frame)
	    self.helix_axis_coords[i,:] = xcoord[:] 
	    self.helix_axis_frames[:,:,i] = xframe[:,:]
        #__for i in xrange(0, num_coords)

        # Reset helix end coordinates.
        self.set_end_coords()

    def get_center(self):
        """ Get helix geometric center. 

            Returns the helix geometric center (NumPy 3x1 array[float]). 
        """
        center = np.mean(self.helix_axis_coords, axis=0) 
        return center 

    def get_domain_ids(self):
        """ Get the list of IDs of the domains in this helix.  """
        domain_ids = set()
        for base in self.staple_bases:
            if base.domain != None:
                domain_ids.add(base.domain)
        for base in self.scaffold_bases: 
            if base.domain != None:
                domain_ids.add(base.domain)
        return list(domain_ids)

    def compute_design_crossovers(self,dna_structure):
        """ Determine the scaffold and staple crossovers in the designed structure.
        """
        self._logger.debug("=================== compute design cross-overs p helix num %d ===================" % self.lattice_num)
        self._logger.debug("Helix polarity %s " % self.scaffold_polarity)
        self._logger.debug("Helix connectivity: %d " % len(self.helix_connectivity)) 
        strands = dna_structure.strands
        for connection in self.helix_connectivity:
            num = connection.to_helix.lattice_num
            self._logger.debug("Crossover helix num %d" % num)
            self._logger.debug("Number of staple bases: %d" % len(self.staple_bases))
            last_crossover = None
            for base_list in [self.staple_bases,self.scaffold_bases]:
                for base in base_list:
                    if not base:
                        continue
                    self._logger.debug(">>> Base  id %d" % base.id)
                    if (base.down != None):
                        if (base.down.h != base.h) and (base.down.h == num):
                            self._logger.debug("base:%4d  p:%4d  h:%4d" % (base.id, base.p, base.h))
                            self._logger.debug("  xd:%4d  p:%4d  h:%4d" % (base.down.id, base.down.p, base.down.h))
                            strand = dna_structure.get_strand(base.strand)
                            crossover = DnaHelixCrossover(self,connection,base,strand)
                            connection.crossovers.append(crossover)
                    #__if (down != -1)

                    if (base.up != None):
                        if (base.up.h != base.h) and (base.up.h == num):
                            self._logger.debug("base:%4d  p:%4d  h:%4d" % (base.id, base.p, base.h))
                            self._logger.debug("  xu:%4d  p:%4d  h:%4d" % (base.up.id, base.up.p, base.up.h))
                            strand = dna_structure.get_strand(base.strand)
                            crossover = DnaHelixCrossover(self,connection,base,strand)
                            connection.crossovers.append(crossover)
                    #__if (up != -1)
                #__for base in base_list
            #__for base_list in [self.staple_bases,self.scaffold_bases]
            self._logger.debug(">>> added %d crossovers " % len(connection.crossovers))
        #__for connection in self.helix_connectivity:
    #__def compute_design_crossovers_p

    def remove_bases(self, base_list):
        """ Remove a list of bases from the helix.

            Arguments:
                base_list (List[DnaBase]): The list of bases to remove.

            This function is used to remove bases after a design file has been processed. 
            For example, it is called when removing strands from a structure.
        """
        self._logger.debug("========== remove bases ==========")
        self._logger.debug("Number of staple bases %d" % len(self.staple_bases))
        self._logger.debug("Number of bases to remove %d" % len(base_list))

        # Create base helix position maps.
        removed_staple_bases = set()
        removed_scaffold_bases = set()
        for base in base_list:
            if base.is_scaf:
                removed_scaffold_bases.add(base.p)
            else:
                removed_staple_bases.add(base.p)
        #__for base in base_list

        # Modify list of staple bases.
        remaining_staple_bases = []
        for base in self.staple_bases:
            if base.p not in removed_staple_bases:
                remaining_staple_bases.append(base)
            else:
                base.up = None
                base.down = None
                if base.across: 
                    base.across.across = None 
                    base.across = None 
            #__if base.p not in removed_staple_bases
        #__for base in self.staple_bases

        self._logger.debug("Number of bases removed %d" % (len(self.staple_bases) - len(remaining_staple_bases)))
        self._logger.debug("Number of bases remaining %d" % len(remaining_staple_bases)) 
        self.staple_bases = remaining_staple_bases
        for base in self.staple_bases:
            self._logger.debug("Base ID %d  h %d  p %d" % (base.id, base.h, base.p))

        # Regenerate the helix coordinate and reference frame  
        # arrays from the modified list of bases.
        self.regenerate_coordinate_arrays()

        # Rebuild base position maps.
        self.build_base_pos_maps()

    def process_base_deletes(self):
        """ Process bases flagged for deletion.

            This function is used to remove bases flagged for deletion (DnaBase.num_deletion) in the 
            design file. It is called when the design file is being processed. Bases are removed from 
            the helix staple_bases and scaffold_bases. The coordinates and frames for the helix are then 
            recreated from the new list of bases.

            A list of deleted bases is returned.
        """
        deleted_bases = []
  
        # Delete staple bases. 
        num_del = 0
        for base in self.staple_bases: 
            if base.num_deletions != 0:
                num_del += 1
        #__for base in self.staple_bases

        if num_del != 0:
            staple_bases = []
            for base in self.staple_bases: 
                if base.num_deletions == 0:
                    staple_bases.append(base)
                else:
                    base.remove()
                    deleted_bases.append(base)
            #__for base in self.staple_bases
            self.staple_bases = staple_bases 
        #__if num_del != 0

        # Delete scaffold bases. 
        num_del = 0
        for base in self.scaffold_bases: 
            if base.num_deletions != 0:
                num_del += 1
        #__for base in self.scaffold_bases

        if num_del != 0:
            scaffold_bases = []
            for base in self.scaffold_bases: 
                if base.num_deletions == 0:
                    scaffold_bases.append(base)
                else:
                    base.remove()
                    deleted_bases.append(base)
            #__for base in self.scaffold_bases
            self.scaffold_bases = scaffold_bases
        #__if num_del != 0

        # Regenerate the helix coordinate and reference frame  arrays 
        # from the new list of bases.
        if len(deleted_bases) != 0:
            self.regenerate_coordinate_arrays()
            self.build_base_pos_maps()
        return deleted_bases 
    #__def process_base_deletes(self)

    def insert_bases(self, insert_bases):
        """ Insert a list of bases into the helix lists of scaffold and staple bases.

            Arguments:
                insert_bases (List[DnaBase]): The list of bases to insert.

            Bases are inserted into the scaffold and staple base lists at the 
            position given in the base (i.e. DnaBase.p).
        """
        self._logger.debug("=================== insert_bases ===================")
        self._logger.debug("Number of bases to insert %d" % len(insert_bases))
        for insert_base in insert_bases:
            if insert_base.is_scaf:
                for i,base in enumerate(self.scaffold_bases): 
                    if base.p == insert_base.p: 
                        self._logger.debug("Insert scaffold base at position %d" % insert_base.p)
                        self._logger.debug("    base coord        %s" % str(base.coordinates))
                        self._logger.debug("    insert base coord %s" % str(insert_base.coordinates))
                        self.scaffold_bases.insert(i,insert_base)
                        break
                #__for i,base in enumerate(self.scaffold_bases)
            else:
                for i,base in enumerate(self.staple_bases): 
                    if base.p == insert_base.p: 
                        self._logger.debug("Insert staple base at position  %d" % insert_base.p)
                        self._logger.debug("    base coord        %s" % str(base.coordinates))
                        self._logger.debug("    insert base coord %s" % str(insert_base.coordinates))
                        self.staple_bases.insert(i,insert_base)
                        break
                #__for i,base in enumerate(self.staple_bases)
        #__for insert_base in insert_bases

        # Regenerate the helix coordinate arrays from the new list of bases.
        if len(insert_bases) != 0:
            self.regenerate_coordinate_arrays()
            self.build_base_pos_maps()
    #__def insert_bases(self, insert_bases)

    def regenerate_coordinate_arrays(self):
        """ Regenerate the coordinate arrays after a change in the base lists.

            This will regenerate the coordinate arrays after deletions or insertions. The base 
            helix position cannot be used to sort the arrays because inserted bases will have the 
            same position. The arrays will therefore be sorted by base distance from the beginning
            helix coordinates.
        """
        self._logger.debug("=================== regenerate_coordinate_arrays %d ===================" % self.lattice_num)

        # Define the origin and axis used to calculate the distance along a helix.
        origin = self.end_coordinates[0]
        point1 = self.end_coordinates[1]
        axis = point1 - origin 

        # Compute the distance for each base from the helix end.
        distances = []
        for base in itertools.chain(self.scaffold_bases,self.staple_bases): 
            coords = base.coordinates;
            v = coords - origin
            d = np.dot(axis,v)
            distances.append((d,base))
        #__for base in itertools.chain(self.scaffold_bases,self.staple_bases)

        # Sort and count the number of unique distances.
        sorted_distances = sorted(distances,key=lambda x: x[0]) 
        num_unique_dist = 0
        last_d = None
        for entry in sorted_distances:
            d = entry[0]
            if d != last_d:
                num_unique_dist += 1
            last_d = d
        #__for entry in sorted_distances

        # Create new coordinates and reference frame arrays.
        axis_coords = np.zeros((num_unique_dist,3), dtype=float)
        axis_frames = np.zeros((3,3,num_unique_dist), dtype=float)
        num_unique_dist = 0
        last_d = None
        for entry in sorted_distances:
            d = entry[0]
            base = entry[1]
            if d != last_d:
                axis_coords[num_unique_dist] = base.coordinates 
                axis_frames[:,:,num_unique_dist] = base.ref_frame
                num_unique_dist += 1
            last_d = d
        #__for entry in sorted_distances
        self.helix_axis_coords = axis_coords
        self.helix_axis_frames = axis_frames
        self.set_end_coords()
    #__def regenerate_coordinate_arrays(self)

#__class DnaStructureHelix


class DnaHelixConnection(object):
    """ This class stores information for a pair of connected helices.

        Attributes:
            from_helix (DnaStructureHelix): The source helix of the connection. 
            to_helix (DnaStructureHelix): The destination helix the connection. 
            direction (NumPy 3x1 adarray[float]): The unit vector in the direction of the connection. 
            crossovers (List[DnaHelixCrossover]): The list of crossovers between the two helices.
    """

    def __init__(self, helix1, helix2):
        self.from_helix = helix1
        self.to_helix = helix2
        self.direction = None
        self._compute_direction()
        self.crossovers = []

    def _compute_direction(self):
        """ Compute the unit vector in the direction of the adjacent helix.
            
            This function uses helix axes to compute the direction so it can be used for lattice or off-lattice geometries.
            However, for lattice-based geometries the directions can be calculated implicitly using a lookup table.
        """
        # Get the first helix axis and a point on that axis from the staple bases. 
        # If there is no staple then use the scaffold.
        helix1 = self.from_helix
        if len(helix1.staple_bases) != 0:
            helix1_base = helix1.staple_bases[0]
        elif len(helix1.scaffold_bases) != 0:
            helix1_base = helix1.scaffold_bases[0]
        pt1 = helix1_base.coordinates
        axis1 = [helix1.end_frames[0,2,0], helix1.end_frames[1,2,0], helix1.end_frames[2,2,0]]

        # Get the second (adjacent) helix axis and a point on that axis.
        helix2 = self.to_helix
        if len(helix2.staple_bases) != 0:
            helix2_base = helix2.staple_bases[0]
        elif len(helix2.scaffold_bases) != 0:
            helix2_base = helix2.scaffold_bases[0]
        pt2 = helix2_base.coordinates
        axis2 = [helix2.end_frames[0,2,0], helix2.end_frames[1,2,0], helix2.end_frames[2,2,0]]
        axis2_length = np.linalg.norm(axis2)

        # Compute the unit vector in the direction of the adjacent helix.
        vec = pt1 - pt2
        d = np.dot(axis2,vec) / axis2_length
        a2pt = pt2 + np.dot(axis2,d)
        self.direction = a2pt - pt1
        self.direction = self.direction / np.linalg.norm(self.direction)
    #__def _compute_direction
#__class DnaHelixConnection


class DnaHelixCrossover(object):
    """ This class stores information for a crossover between two helices.

        Attributes:
            helix (DnaStructureHelix): The helix the crossover is in.
            helix_connection (DnaHelixConnection): The connection object defining the crossover helix.
            crossover_base (Base): The base where the crossover occurs. 
            strand (DnaStrand): The strand the crossover is in.
    """
    def __init__(self, helix, helix_connection, crossover_base, strand):
        self.helix = helix
        self.helix_connection = helix_connection
        self.crossover_base = crossover_base
        self.strand = strand

#__class DnaHelixCrossover
