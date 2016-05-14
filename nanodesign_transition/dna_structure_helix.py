#!/usr/bin/env python
"""
This module defines the classes used to define a structure helix of a DNA structure.

A structure helix is a region in a DNA structure that forms a cylindrical structural element. It can be composed of 
one or two DNA strands. 

"""
from collections import OrderedDict
import inspect
import json
import logging
import numpy as np
from parameters import DnaParameters
from converters.cadnano.common import CadnanoLatticeType

# temp code to handle objects as they are being transitioned into the main package
try:
    # TODO: JS 3/25 This will need to change at some point once everything is transitioned.
    import os.path
    import sys
    base_path = os.path.abspath( os.path.dirname(__file__) + '/../' )
    sys.path.append(base_path)
    import nanodesign as nd
    from nanodesign_transition.lattice import Lattice
    sys.path = sys.path[:-1]
except ImportError:
    print "Cannot locate nanodesign package, it hasn't been installed in main packages, and is not reachable relative to the nanodesign_transition directory."
    raise ImportError

class DnaStructureHelix(object):
    """ This class stores information for a DNA structure helix. 

        Attributes:
            end_coordinates (NumPy 2x3 ndarray[float]): The coordinates at the ends of the helix.
            end_frames (NumPy 3x3xN ndarray[float]): The coordinate frames at the ends of the helix.
            helix_axis_frames (NumPy 3x3xN ndarray[float]): The coordinate frames of base nodes along the helix axis, 
                where N is the size of the caDNAno virtual helix.
            helix_axis_nodes (NumPy Nx3 ndarray[float]): The coordinates of base nodes along the helix axis, where N is the 
                size of the caDNAno virtual helix.
            helix_connectivity (List[DnaHelixConnection]): The list helices that this helix is connected to.
            id (int): Helix ID (1 - number of helices in a structure).
            lattice_col (int): The caDNAno lattice column number.
            lattice_num (int): The caDNAno helix number.
            lattice_row (int): The caDNAno lattice row number. 
            possible_scaffold_crossovers (List[(DnaStructureHelix,int)]): The list of possible scaffold crossovers for this 
                helix. Each entry in the list is a 2-tuple giving the helix it crosses over to (object pointer) and the helix 
                position where the crossover occurs.
            possible_staple_crossovers (List[(DnaStructureHelix,int)]): The list of possible staple crossovers for this 
                helix. Each entry in the list is a 2-tuple giving the helix it crosses over to (object pointer) and the 
                helix position where the crossover occurs.
            scaffold_base_list (List[DnaBase]): The list storing helix scaffold bases. This list is the same size for all helices
                so it may have None entries where scaffold strands do not have a base at that location. 
            scaffold_polarity (DnaPolarity): The polarity of the scaffold. 
            staple_base_list (List[DnaBase]): The list storing helix staple bases. This list is the same size for all helices
                so it may have None entries where staple strandis do not have a base at that location. 
            start_strand (int): The start position of the scaffold or staple strand in the helix. This used to index into the
                helix_axis_coordinates and helix_axis_frames arrays.  
            end_strand (int): The end position of the scaffold or staple strand in the helix. This used to index into the
                helix_axis_coordinates and helix_axis_frames arrays.  
    """ 

    # TODO (DaveP) We need to remove the references to caDNAno lattice-based information.

    def __init__(self, id, scaffold_polarity, helix_axis_nodes, helix_axis_frames, start_strand, end_strand):
        """ Initialize a DnaStructureHelix object.

            Arguments:
                id (int): The helix ID.
                scaffold_polarity (DnaPolarity): The helix polarity. 
                helix_axis_nodes (NumPy Nx3 ndarray[float]): The coordinates of base nodes along the helix axis. 
                helix_axis_frames (NumPy 3x3xN ndarray[float]): The coordinate frames base nodes alonge he helix axis.
                start_strand (int): The start position of the scaffold or staple strand in the helix. 
                end_strand (int): The end position of the scaffold or staple strand in the helix. 
        """
        self.id = id
        self.helix_axis_frames = helix_axis_frames
        self.helix_axis_nodes = helix_axis_nodes
        self.scaffold_polarity = scaffold_polarity 
        self.start_strand = start_strand 
        self.end_strand = end_strand
        self.lattice_row = -1
        self.lattice_col = -1
        self.lattice_num = -1
        self.staple_base_list = []
        self.scaffold_base_list = []
        self.helix_connectivity = []
        self.possible_staple_crossovers = []
        self.possible_scaffold_crossovers = []

        # Set helix ends coordinates.
        self.end_coordinates = np.zeros((2,3), dtype=float)
        self.end_frames = np.zeros((3,3,2), dtype=float)
        if self.start_strand != -1:
            self.end_coordinates[0] = self.helix_axis_nodes[start_strand]
            self.end_coordinates[1] = self.helix_axis_nodes[end_strand]
            self.end_frames[:,:,0] = helix_axis_frames[:,:,start_strand]
            self.end_frames[:,:,1] = helix_axis_frames[:,:,end_strand]
        else:
            self.end_coordinates[0] = None
            self.end_coordinates[1] = None
            self.end_frames[:,:,0] = None
            self.end_frames[:,:,1] = None

    def _setup_logging(self,name):
        """ Set up logging."""
        logger = logging.getLogger(name)
        logger.setLevel(logging.INFO)
        # create console handler and set format
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        return logger

    def get_start_pos(self):
        """ Get the starting helix position of the scaffold or staple strands. 
            This is used for visualization.
        """
        num_bases = len(self.staple_base_list)
        staple_start_pos = next((base.p for base in self.staple_base_list if base != None),num_bases)
        scaffold_start_pos = next((base.p for base in self.scaffold_base_list if base != None),num_bases)
        start_pos = min(staple_start_pos, scaffold_start_pos)
        return start_pos 

    def get_domain_ids(self):
        """ Get the list of IDs of the domains in this helix.  """
        domain_ids = set()
        for base in self.staple_base_list:
            if base:
                domain_ids.add(base.domain)
        for base in self.scaffold_base_list: 
            if base:
                domain_ids.add(base.domain)
        return list(domain_ids)

    def compute_design_crossovers(self,dna_structure):
        """ Determine the scaffold and staple crossovers in the designed structure.
        """
        logger = self._setup_logging(__name__+str(self.id))
        #logger.setLevel(logging.DEBUG)
        logger.debug("=================== compute design cross-overs helix num %d ===================" % self.lattice_num)
        logger.debug(">>> helix polarity %s " % self.scaffold_polarity)
        logger.debug(">>> helix connectivity: %d " % len(self.helix_connectivity)) 
        base_connectivity = dna_structure.base_connectivity
        strands = dna_structure.strands
        for connection in self.helix_connectivity:
            num = connection.to_helix.lattice_num
            logger.debug(">>> crossover helix num %d" % num)
            logger.debug(">>> num staple bases: %d" % len(self.staple_base_list))
            last_crossover = None
            for base_list in [self.staple_base_list,self.scaffold_base_list]:
                for base in base_list:
                    if not base:
                        continue
                    down = base.down
                    up = base.up

                    if (down != -1):
                        down_base = base_connectivity[down-1]
                        if (down_base.h != base.h) and (down_base.h == num):
                            logger.debug("base:%4d  p:%4d  h:%4d" % (base.id, base.p, base.h))
                            logger.debug("  xd:%4d  p:%4d  h:%4d" % (down_base.id, down_base.p, down_base.h))
                            strand = dna_structure.get_strand(base.strand)
                            crossover = DnaHelixCrossover(self,connection,base,strand)
                            connection.crossovers.append(crossover)
                    #__if (down != -1)

                    if (up != -1):
                        up_base = base_connectivity[up-1]
                        if (up_base.h != base.h) and (up_base.h == num):
                            logger.debug("base:%4d  p:%4d  h:%4d" % (base.id, base.p, base.h))
                            logger.debug("  xu:%4d  p:%4d  h:%4d" % (up_base.id, up_base.p, up_base.h))
                            strand = dna_structure.get_strand(base.strand)
                            crossover = DnaHelixCrossover(self,connection,base,strand)
                            connection.crossovers.append(crossover)
                    #__if (up != -1)

                #__for base in self.base_list
            #__for base_list in [self.staple_base_list,self.scaffold_base_list]
            logger.debug(">>> added %d crossovers " % len(connection.crossovers))
        #__for connection in self.helix_connectivity:
    #__def compute_design_crossovers
#__class DnaStructureHelix


class DnaHelixConnection(object):
    """ This class stores information for a pair of connected helices.

        Attributes:
            from_helix (DnaStructureHelix): The source helix of the connection. 
            to_helix (DnaStructureHelix): The destination helix the connection. 
            direction (NumPy 3x1 ndarray[float]): The unit vector in the direction of the connection. 
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
        start_pos = next((i for i in xrange(0,len(helix1.staple_base_list)) if helix1.staple_base_list[i] != None),-1)
        if start_pos == -1: 
            start_pos = next((i for i in xrange(0,len(helix1.scaffold_base_list)) if helix1.scaffold_base_list[i] != None),-1)
            helix1_base = helix1.scaffold_base_list[start_pos]
        else:
            helix1_base = helix1.staple_base_list[start_pos]
        pt1 = helix1.helix_axis_nodes[helix1_base.p]
        axis1 = [helix1.end_frames[0,2,0], helix1.end_frames[1,2,0], helix1.end_frames[2,2,0]]

        # Get the second (adjacent) helix axis and a point on that axis.
        helix2 = self.to_helix
        start_pos = next((i for i in xrange(0,len(helix2.staple_base_list)) if helix2.staple_base_list[i] != None),-1)
        if start_pos == -1: 
            start_pos = next((i for i in xrange(0,len(helix2.scaffold_base_list)) if helix2.scaffold_base_list[i] != None),-1)
            helix2_base = helix2.scaffold_base_list[start_pos]
        else:
            helix2_base = helix2.staple_base_list[start_pos]
        #print("    start pos %d" % start_pos) 
        pt2 = helix2.helix_axis_nodes[helix2_base.p]
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
