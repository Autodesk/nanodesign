#!/usr/bin/env python
"""
This module defines the classes used to define the connectivity and geometry of a DNA structure. 

A DNA structure consists of a number of scaffold and staple strands (DNA origami), or oligo strands alone, bound together 
to form a designed geometric shape.
"""
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
    sys.path = sys.path[:-1]
except ImportError:
    print "Cannot locate nanodesign package, it hasn't been installed in main packages, and is not reachable relative to the nanodesign_transition directory."
    raise ImportError


class DnaStructure(object):
    """ This class stores the base connectivity and geometry for a DNA model. 

        Attributes:
            base_connectivity (DnaBase): A list of DnaBase objects.
            helix_axis_nodes (numpy Nx3 float darray): The coordinates of base nodes along a helix axis. 
            helix_axis_frames (numpy 3x3xN float darray): The reference frames of bases along a helix axis. The reference frame
                    is a right-handed coordinate frame (e1,e2,e3) attached to each base. e1 points in the direction of the major 
                    groove, e2 runs along the long helix axis and e3 = e1 x e2.
            strands (DnaStrand): A list a DnaStrand objects. 
    """ 

    def __init__(self, name="dna structure"):
        self.name = name
        self.base_connectivity = None
        self.helix_axis_nodes = None
        self.helix_axis_frames = None
        self.strands = None
        self.id_nt = None
        self.structure_helices = []
        self.structure_helices_map = dict()
        self.parameters = DnaParameters()
        self.staple_colors = []
        self.domain_list = []
        self.dnodes_map = dict()
        self.lattice_type = CadnanoLatticeType.none
        self._logger = self._setup_logging()

    def _setup_logging(self):
        """ Set up logging."""
        logger = logging.getLogger('dna_structure')
        logger.setLevel(logging.INFO)
        # create console handler and set format
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        return logger

    def add_structure_helices(self, structure_helices):
        """ Add a list of structural helices. """
        for helix in structure_helices:
            self.structure_helices.append(helix)
            self.structure_helices_map[helix.lattice_num] = helix

    def compute_aux_data(self):
        """ Compute auxiallry data. """
        for strand in self.strands:
            strand.dna_structure = self
        self._set_strand_helix_references()
        self._set_helix_bases()
        self._compute_strand_helix_references()
        self._compute_domains()

    def get_domains(self):
        if (not self.domain_list): 
            self._compute_domains()
        return self.domain_list

    def get_strand(self,id):
        """ Get a strand from an id. """
        for strand in self.strands:
            if (strand.id == id):
               return strand
        #__for strand in self.strands
        return None




    def _set_helix_bases(self):
        """ Set the bases for a helix. """
        for helix in self.structure_helices:
            num = helix.lattice_num
            hsize = len(helix.helix_axis_nodes)
            staple_base_list = [None]*hsize
            scaffold_base_list = [None]*hsize
            for base in self.base_connectivity:
                if (base.h != num):
                    continue
                if base.is_scaf:
                    scaffold_base_list[int(base.p)] = base
                else:
                    staple_base_list[int(base.p)] = base
            helix.staple_base_list = staple_base_list
            helix.scaffold_base_list = scaffold_base_list

    def _set_strand_helix_references(self):
        """ Set the helices referenced by each strand. """
        for strand in self.strands:
            for id in strand.tour:
                base = self.base_connectivity[id-1]
                helix = self.structure_helices_map[base.h]
                strand.add_helix(helix)
        #__for strand in self.strands__

    def _compute_domains(self):
        """ Compute DNA domains. 
      
            Domains are computed using the list of bases defined for each structural helix in the model using the base
            position within the helix. Helix bases are traversed from index 0 regardless of the helix polarity. The list 
            storing helix bases (staple and scaffold) are the same size for all helices and contain None where bases are 
            not defined. 

            The boolean array strand_breaks[] is used to mark the location (relative position in the helix) of where a 
            strand changes its binding to another helix (cross-over), becomes a single strand or begins/ends. Entries in 
            strand_breaks[] are set to True for each base position in staple and scaffold strands within a helix that 
            represent a binding change. Pairs of entries in strand_breaks[] are then used to define the boundaries of domains.
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("===================== compute domains =====================")
        num_domains = 0
        self.domain_list = []

        for helix in self.structure_helices:
            num = helix.lattice_num
            self._logger.debug("---------- helix num %d ---------- " % num)
            staple_base_list = helix.staple_base_list 
            scaffold_base_list = helix.scaffold_base_list
            helix_size = len(scaffold_base_list)
            strand_breaks = [False]*helix_size
            if helix.scaffold_polarity == "3'": 
                three_to_five = True 
            else:
                three_to_five = False

            # Add cross-overs for staples.
            strand_breaks = self._set_strand_breaks(strand_breaks, staple_base_list)

            # Add cross-overs for scaffolds.
            strand_breaks = self._set_strand_breaks(strand_breaks, scaffold_base_list)

            # Create domains from pairs of True entries in strand_breaks[].
            self._logger.debug(">>> create domains:")
            scaffold = False

            # First iterate over staple bases, then scaffold bases.
            for base_list in [staple_base_list, scaffold_base_list]:
                start = -1
                # Set how bases are added to a domain: in order or reverse order.
                if (scaffold):
                    reverse = three_to_five 
                else:
                    reverse = not three_to_five

                for i in xrange(0,len(strand_breaks)):
                    base = base_list[i]
                    if (not (strand_breaks[i] and base)):
                        continue 
                    if (start != -1):
                        self._logger.debug("    domain: %4d  start pos:%4d  end pos: %4d "  % (num_domains,start,i))
                        domain = nd.Domain(num_domains,helix)
                        strand = self.get_strand(base.strand)
                        domain.strand = strand
                        domain.color = strand.color
                        helix.domain_list.append(domain)
                        self.domain_list.append(domain)
                        num_domains += 1
                        for j in xrange(start,i+1):
                            base = base_list[j]
                            if (reverse):
                                domain.base_list.insert(0,base)
                            else:
                                domain.base_list.append(base)
                            base.domain = domain.id
                        start = -1
                    else:
                        start = i
                #__for i in xrange(0,len(strand_breaks)):
                scaffold = True
            #__for base_list in [staple_base_list, scaffold_base_list]
        self._logger.debug(">>> Created %d domains." % num_domains)

        # set the strand and domain each domain is connected to.
        for domain in self.domain_list:
            across = -1
            for base in domain.base_list:
                if (base.across != -1):
                    across = base.across
                    break
            #__for base in domain.base_list
            conn_dom = -1
            conn_strand = -1
            if (across != -1):
                across_base = self.base_connectivity[across-1]
                conn_dom = across_base.domain
                conn_strand = across_base.strand
            domain.connected_strand = conn_strand
            domain.connected_domain = conn_dom
    #__def _compute_domains(self)

    def _set_strand_breaks(self, strand_breaks, base_list):
        """ Set the location of a change in strand binding.

            A strand break occurs if a strand:
                1) crosses over to another helix 
                2) ends/begins                   
                3) becomes a single strand      

            Arguments:
                strand_breaks (list[bool]): The locations of strand binding change. 
                base_list (list[DnaBase]): The list bases in a helix. 

            Returns:
                strand_breaks (list[bool]): The locations of strand binding change. 
        """
        # Find the location of the first base in the helix.
        start_pos = next((i for i in xrange(0,len(base_list)) if base_list[i] != None),-1)
        base = base_list[start_pos]
        across_base = base.across
        curr_across_base_sign = np.sign(base.across)
        strand_breaks[start_pos] = True

        # Set the locations of binding changes for each base.
        for i in xrange(start_pos+1,len(base_list)):
            base = base_list[i]
            if (not base):
                continue
            across_base_sign = np.sign(base.across)

            # Check for a cross-over or end/start of a strand.
            if self._check_base_crossover(base):
                curr_across_base_sign = np.sign(base.across)
                strand_breaks[i] = True

            # A change in the base across sign signals that a single strand is starting/ending.
            elif (across_base_sign != curr_across_base_sign):
                 curr_across_base_sign = np.sign(base.across)
                 strand_breaks[i-1] = True
                 strand_breaks[i] = True
        #__for i in xrange(start_pos+1,len(base_list)):
        return strand_breaks 

    def _check_base_crossover(self,base):
        """ Check if a base marks a cross-over to another helix or the end/start of a strand.

            A -1 for a bases's down/up signifies that a strand is starting or ending. A change in 
            a base's up/down helix signifies a cross-over. 

            Arguments:
                base (DnaBase): The base to check for a cross-over. 

            Returns: True if the base crosses over to another helix or ends/begins.
        """
        down = base.down
        if (down == -1):
            return True

        down_base = self.base_connectivity[down-1]
        if (down_base.h != base.h):
            return True

        up = base.up
        if (up == -1):
            return True
        up_base = self.base_connectivity[up-1]

        if (up_base.h != base.h):
            return True

        return False

    def _compute_strand_helix_references(self):
        """ Set the virtual helices referenced by each strand. """
        for strand in self.strands:
            for id in strand.tour:
                base = self.base_connectivity[id-1]
                helix = self.structure_helices_map[base.h]
                strand.add_helix(helix)
        #__for strand in self.strands__


#__class DnaStructure(object):


class DnaStructureHelix(object):
    """ This class stores information for a DNA structure helix. 

        A structure helix is a region in a DNA structure that forms a cylindrical structural element. It can be composed of 
        one or two DNA strands. 

        Attributes:
            id (int): Helix ID (1 - number of helices in a structure).
            staple_base_list (list[DnaBase]): The list storing helix staple bases. This list is the same size for all helices. 
            scaffold_base_list (list[DnaBase]): The list storing helix scaffold bases. This list is the same size for all helices. 
    """ 
    def __init__(self, id):
        self.id = id
        self.lattice_row = -1
        self.lattice_col = -1
        self.lattice_num = -1
        self.staple_colors = []
        self.end_coordinates = np.zeros((2,3), dtype=float)
        self.end_frames = np.zeros((3,3,2), dtype=float)
        self.helix_axis_nodes = None
        self.domain_list = []
        self.scaffold_polarity = "3'"
        self.staple_base_list = []
        self.scaffold_base_list = []

        # note [davep] are we using these?
        self.start_pos = -1
        self.triad = None
        self.id_nt = None
        self.helix_topology = None



