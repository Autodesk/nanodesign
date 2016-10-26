#!/usr/bin/env python
"""
This module defines the classes used to define a structure helix of a DNA structure.

A structure helix is a region in a DNA structure that forms a cylindrical structural element. It can be composed of 
one or two DNA strands. 

"""
from collections import OrderedDict
import inspect
import itertools
import json
import logging
import numpy as np

# Within package imports
from .parameters import DnaParameters
from ..converters.cadnano.common import CadnanoLatticeType
from .lattice import Lattice


# # temp code to handle objects as they are being transitioned into the main package
# try:
#     # TODO: JS 3/25 This will need to change at some point once everything is transitioned.
#     import os.path
#     import sys
#     base_path = os.path.abspath( os.path.dirname(__file__) + '/../' )
#     sys.path.append(base_path)
#     import nanodesign as nd
#     from nanodesign_transition.lattice import Lattice
#     sys.path = sys.path[:-1]
# except ImportError:
#     print "Cannot locate nanodesign package, it hasn't been installed in main packages, and is not reachable relative to the nanodesign_transition directory."
#     raise ImportError

class DnaStructureHelix(object):
    """ This class stores information for a DNA structure helix. 

        Attributes:
            count (int): The helix count (0 based) of the order the helix was processed. This is needed when
                writing a caDNAno file.
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
            possible_scaffold_crossovers (List[(DnaStructureHelix,int)]): The list of possible scaffold crossovers for this 
                helix. Each entry in the list is a 2-tuple giving the helix it crosses over to (object pointer) and the helix 
                position where the crossover occurs.
            possible_staple_crossovers (List[(DnaStructureHelix,int)]): The list of possible staple crossovers for this 
                helix. Each entry in the list is a 2-tuple giving the helix it crosses over to (object pointer) and the 
                helix position where the crossover occurs.
            scaffold_bases (List[DnaBase]): The list of helix scaffold bases. 
            scaffold_polarity (DnaPolarity): The polarity of the scaffold. 
            staple_bases (List[DnaBase]): The list of helix staple bases. 

        The geometry of a helix is determined by the coordinates of base nodes along its axis. The coordinates and reference
        frames of the bases defined for the helix are references into the helix_axis_coords and helix_axis_frames arrays.
    """ 

    # TODO (DaveP) We need to remove the references to caDNAno lattice-based information.

    def __init__(self, count, id, scaffold_polarity, helix_axis_coords, helix_axis_frames, scaffold_bases, staple_bases):
        """ Initialize a DnaStructureHelix object.

            Arguments:
                count (int): The helix count (0 based) of the order the helix was processed.
                id (int): The helix ID.
                scaffold_polarity (DnaPolarity): The helix polarity. 
                helix_axis_coords (NumPy Nx3 ndarray[float]): The coordinates of base nodes along the helix axis. 
                helix_axis_frames (NumPy 3x3xN ndarray[float]): The coordinate frames base nodes alonge the helix axis.
                scaffold_bases (List[DnaBase]): The list of scaffold bases defined for the helix.
                staple_bases (List[DnaBase]): The list of staple bases defined for the helix.
        """
        self.id = id
        self.count = count 
        self.staple_bases = staple_bases
        self.scaffold_bases = scaffold_bases
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
        self.logger = self._setup_logging()
        # Set helix ends coordinates.
        self.end_coordinates = np.zeros((2,3), dtype=float)
        self.end_frames = np.zeros((3,3,2), dtype=float)
        self.set_end_coords()

    def _setup_logging(self):
        """ Set up logging. """
        logger = logging.getLogger(__name__ + ":" + str(self.id))
        logger.setLevel(logging.INFO)
        # Create console handler and set format.
        if not len(logger.handlers):
            console_handler = logging.StreamHandler()
            formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
            console_handler.setFormatter(formatter)
            logger.addHandler(console_handler)
        return logger

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

    def get_start_pos(self):
        """ Get the starting helix position of the scaffold or staple strands. 
            This is used for visualization.
        """
        num_bases = len(self.staple_bases)
        staple_start_pos = next((base.p for base in self.staple_bases if base != None),num_bases)
        scaffold_start_pos = next((base.p for base in self.scaffold_bases if base != None),num_bases)
        start_pos = min(staple_start_pos, scaffold_start_pos)
        return start_pos 

    def apply_xform(self, xform):
        """ Apply a transformation to the helix coordiates and reference frames.

            Arguments:
                xform (Xform): The transformation to apply to the helx geometry.
        """
        self.logger.setLevel(logging.INFO)
        #self.logger.setLevel(logging.DEBUG)
        self.logger.debug("=================== apply xform to helix %d ===================" % self.id)
        R = xform.rotation_matrix
        #xform.print_transformation()
        translation = xform.translation 
        center = xform.center 
        coord = np.array([0.0,0.0,0.0], dtype=float) 
        frame = np.zeros((3,3), dtype=float)
        xcoord = np.array([0.0,0.0,0.0], dtype=float) 
        xframe = np.zeros((3,3), dtype=float)
        num_coords = len(self.helix_axis_coords)
        self.logger.debug("Number of coordinates %d" % num_coords) 
        self.logger.debug("Xform center (%g %g %g)" % (center[0], center[1], center[2])) 
        self.logger.debug("Xform translation (%g %g %g)" % (translation[0], translation[1], translation[2])) 
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
        self.logger.setLevel(logging.INFO)
        #self.logger.setLevel(logging.DEBUG)
        self.logger.debug("=================== compute design cross-overs p helix num %d ===================" % self.lattice_num)
        self.logger.debug("Helix polarity %s " % self.scaffold_polarity)
        self.logger.debug("Helix connectivity: %d " % len(self.helix_connectivity)) 
        strands = dna_structure.strands
        for connection in self.helix_connectivity:
            num = connection.to_helix.lattice_num
            self.logger.debug("Crossover helix num %d" % num)
            self.logger.debug("Number of staple bases: %d" % len(self.staple_bases))
            last_crossover = None
            for base_list in [self.staple_bases,self.scaffold_bases]:
                for base in base_list:
                    if not base:
                        continue
                    self.logger.debug(">>> Base  id %d" % base.id)
                    if (base.down != None):
                        if (base.down.h != base.h) and (base.down.h == num):
                            self.logger.debug("base:%4d  p:%4d  h:%4d" % (base.id, base.p, base.h))
                            self.logger.debug("  xd:%4d  p:%4d  h:%4d" % (base.down.id, base.down.p, base.down.h))
                            strand = dna_structure.get_strand(base.strand)
                            crossover = DnaHelixCrossover(self,connection,base,strand)
                            connection.crossovers.append(crossover)
                    #__if (down != -1)

                    if (base.up != None):
                        if (base.up.h != base.h) and (base.up.h == num):
                            self.logger.debug("base:%4d  p:%4d  h:%4d" % (base.id, base.p, base.h))
                            self.logger.debug("  xu:%4d  p:%4d  h:%4d" % (base.up.id, base.up.p, base.up.h))
                            strand = dna_structure.get_strand(base.strand)
                            crossover = DnaHelixCrossover(self,connection,base,strand)
                            connection.crossovers.append(crossover)
                    #__if (up != -1)
                #__for base in base_list
            #__for base_list in [self.staple_bases,self.scaffold_bases]
            self.logger.debug(">>> added %d crossovers " % len(connection.crossovers))
        #__for connection in self.helix_connectivity:
    #__def compute_design_crossovers_p

    def remove_bases(self, base_list):
        """ Remove a list of bases from the helix.

            Arguments:
                base_list (List[DnaBase]): The list of bases to remove.
        """
        self.logger.setLevel(logging.DEBUG)
        self.logger.debug("========== remove bases ==========")
        self.logger.debug("Number of staple bases %d" % len(self.staple_bases))
        self.logger.debug("Number of bases to remove %d" % len(base_list))

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
                #base.remove()
                if base.across: 
                    base.across.across = None 
            #__if base.p not in removed_staple_bases
        #__for base in self.staple_bases

        self.logger.debug("Number of bases removed %d" % (len(self.staple_bases) - len(remaining_staple_bases)))
        self.logger.debug("Number of bases remaining %d" % len(remaining_staple_bases)) 
        self.staple_bases = remaining_staple_bases
        for base in self.staple_bases:
            self.logger.debug("Base ID %d  h %d  p %d" % (base.id, base.h, base.p))

        # Regenerate the helix coordinate and reference frame  
        # arrays from the modified list of bases.
        self.regenerate_coordinate_arrays()

    def process_base_deletes(self):
        """ Process bases flagged for deletion.

            Bases flagged for deletion (Base.skip = 0) are removed from the helix 
            staple_bases and scaffold_bases. The coordinates and frames for the
            helix are then recreated from the new list of bases.

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
        return deleted_bases 
    #__def process_base_deletes(self)

    def insert_bases(self, insert_bases):
        """ Insert a list of bases into the helix lists of scaffold and staple bases.

            Arguments:
                insert_bases (List[DnaBase]): The list of bases to insert.

            Bases are inserted into the scaffold and staple base lists at the 
            position given in the base (i.e. DnaBase.p).
        """
        #self.logger.setLevel(logging.DEBUG)
        self.logger.setLevel(logging.INFO)
        self.logger.debug("=================== insert_bases ===================")
        self.logger.debug("Number of bases to insert %d" % len(insert_bases))
        for insert_base in insert_bases:
            if insert_base.is_scaf:
                for i,base in enumerate(self.scaffold_bases): 
                    if base.p == insert_base.p: 
                        self.logger.debug("Insert scaffold base at position %d" % insert_base.p)
                        self.logger.debug("    base coord        %s" % str(base.coordinates))
                        self.logger.debug("    insert base coord %s" % str(insert_base.coordinates))
                        self.scaffold_bases.insert(i,insert_base)
                        break
                #__for i,base in enumerate(self.scaffold_bases)
            else:
                for i,base in enumerate(self.staple_bases): 
                    if base.p == insert_base.p: 
                        self.logger.debug("Insert staple base at position  %d" % insert_base.p)
                        self.logger.debug("    base coord        %s" % str(base.coordinates))
                        self.logger.debug("    insert base coord %s" % str(insert_base.coordinates))
                        self.staple_bases.insert(i,insert_base)
                        break
                #__for i,base in enumerate(self.staple_bases)
        #__for insert_base in insert_bases

        # Regenerate the helix coordinate arrays from the new list of bases.
        if len(insert_bases) != 0:
            self.regenerate_coordinate_arrays()
        self.logger.setLevel(logging.INFO)
    #__def insert_bases(self, insert_bases)

    def regenerate_coordinate_arrays(self):
        """ Regenerate the coordinate arrays after a change in the base lists.

            This will regenerate the coordinate arrays after deletions or insertions. The base 
            helix position cannot be used to sort the arrays because inserted bases will have the 
            same position. The arrays will therefore be sorted by base distance from the beginning
            helix coordinates.
        """
        #self.logger.setLevel(logging.DEBUG)
        self.logger.setLevel(logging.INFO)
        self.logger.debug("=================== regenerate_coordinate_arrays %d ===================" % self.lattice_num)

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
        self.logger.setLevel(logging.INFO)
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
        else:
            helix1_base = helix1.staple_bases[0]
        pt1 = helix1_base.coordinates
        axis1 = [helix1.end_frames[0,2,0], helix1.end_frames[1,2,0], helix1.end_frames[2,2,0]]

        # Get the second (adjacent) helix axis and a point on that axis.
        helix2 = self.to_helix
        if len(helix2.staple_bases) != 0:
            helix2_base = helix2.staple_bases[0]
        else:
            helix2_base = helix2.staple_bases[0]
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
