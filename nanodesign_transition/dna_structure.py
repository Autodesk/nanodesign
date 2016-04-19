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
        self.strands_map = dict()
        self.id_nt = None
        self.structure_helices = []
        self.structure_helices_map = dict()
        self.structure_helices_coord_map = dict()
        self.parameters = DnaParameters()
        self.staple_colors = []
        self.domain_list = []
        self.dnodes_map = dict()
        self.lattice_type = CadnanoLatticeType.none
        self.lattice = None
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
            self.structure_helices_coord_map[(helix.lattice_row,helix.lattice_col)] = helix

    def compute_aux_data(self):
        """ Compute auxiallry data. """
        for strand in self.strands:
            strand.dna_structure = self
        self._set_strand_helix_references()
        self._set_helix_bases()
        self._compute_strand_helix_references()
        self._compute_domains()
        self._set_helix_connectivity()
        self._compute_helix_design_crossovers()

    def get_domains(self):
        if (not self.domain_list): 
            self._compute_domains()
        return self.domain_list

    def get_strand(self,id):
        """ Get a strand from an id. """
        if not self.strands_map:
            for strand in self.strands:
                self.strands_map[strand.id] = strand 
        if id not in self.strands_map:
            self._logger.error("Failed to find strand id %d." % id)
            return None
        return self.strands_map[id]

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
                    scaffold_base_list[base.p] = base
                else:
                    staple_base_list[base.p] = base
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
        """ Compute DNA domains from strands. 

            Domains are computed by traversing the bases in the scaffold and staple strands of a structure. 
            Domain are bounded by a single->double or double->single strand transitions, crossovers between 
            helices or strand termination. 

            Domains are created using an integer ID starting from 0. 
            Domain objects are stored in self.domain_list[]. A list of domains is also created for each strand.
        """ 
        self._logger.setLevel(logging.INFO)
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("===================== compute domains =====================")
        domain_id = 0
        self.domain_list = []

        # Iterate over the scaffold and staple strands of a structure. 
        for strand in self.strands:
            self._logger.debug("")
            if ( strand.is_scaffold):
                self._logger.debug("==================== scaffold strand %d ====================" % strand.id)
            else:
                self._logger.debug("==================== staple strand %d ====================" % strand.id)

            start_base = self.base_connectivity[strand.tour[0]-1]
            end_base = self.base_connectivity[strand.tour[-1]-1]
            self._logger.debug("Strand number of bases: %3d" % len(strand.tour))
            self._logger.debug("Strand start: h: %3d  p: %3d" % (start_base.h, start_base.p))
            self._logger.debug("Strand end: h: %3d  p: %3d" % (end_base.h, end_base.p))

            # Initialize the domain base list.
            id = strand.tour[0]
            base = self.base_connectivity[id-1]
            curr_across_sign = np.sign(base.across)
            domain_bases = [ base ]

            # Traverse the bases in a strand and create domains.
            for i in xrange(1,len(strand.tour)):
                id = strand.tour[i]
                base = self.base_connectivity[id-1]
                across_sign = np.sign(base.across)
                add_curr_base = True

                # If no domain bases then just continue after checking for sign change.
                if len(domain_bases) == 0:
                    domain_bases.append(base)
                    if curr_across_sign != across_sign:
                        curr_across_sign = across_sign
                    continue

                # Check for a single->double or double->single strand transition.
                if curr_across_sign != across_sign:
                    self._add_domain(domain_id, strand, domain_bases, "sign change")
                    domain_id += 1
                    domain_bases = []
                    curr_across_sign = across_sign

                # Check for a crossover between helices for this base.
                elif self._check_base_crossover(base):
                    last_base = domain_bases[-1]
                    # Make sure the current base is in the same helix.
                    if base.h == last_base.h:
                        domain_bases.append(base)
                        add_curr_base = False
                    self._add_domain(domain_id, strand, domain_bases, "base crossover")
                    domain_id += 1
                    domain_bases = []

                # Check the base paired to this base for: crossover or termination.
                elif base.across > 0:
                    abase = self.base_connectivity[base.across-1]
                    if self._check_base_crossover(abase):
                        domain_bases.append(base)
                        add_curr_base = False
                        self._add_domain(domain_id, strand, domain_bases, "abase crossover")
                        domain_id += 1
                        domain_bases = []
                    # If a strand terminates make sure the current base is in the same strand. 
                    elif (abase.down == -1) or (abase.up == -1):
                        last_base = domain_bases[-1]
                        if last_base.across != -1:
                            last_abase = self.base_connectivity[last_base.across-1]
                            if abase.strand == last_abase.strand:
                                domain_bases.append(base)
                                add_curr_base = False
                        else:
                            domain_bases.append(base)
                            add_curr_base = False
                        self._add_domain(domain_id, strand, domain_bases, "abase start/end")
                        domain_id += 1
                        domain_bases = []
                #__if curr_across_sign != across_sign

                # Add the current base to the current list of domain bases.
                if add_curr_base:
                    domain_bases.append(base)
            #__for i in xrange(1,len(strand.tour))

            # Add a domain for any remaining bases.
            if len(domain_bases) != 0:
                self._add_domain(domain_id, strand, domain_bases, "remaining")
                domain_id += 1
        #__for strand in self.strands__

        self._logger.info("Number of domains computed: %d " % len(self.domain_list))

        # Check if the computed domains are consistent with the strands they were computed from.
        self.check_domains()

        # Set the strand and domain each domain is connected to.
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
        #__for domain in self.domain_list

    def _check_base_crossover(self, base):
        """ Check if there is a crossover to a different helix at the given base.

            Arguments:
                base (DnaBase): The base to check for a crossover.

            Returns True if there is a crossover at the base.
        """
        down = base.down
        if base.down == -1:
            return False
        if self.base_connectivity[base.down-1].h != base.h:
            return True

        if base.up == -1:
            return False
        if self.base_connectivity[base.up-1].h != base.h:
            return True

        return False

    def check_domains(self):
        """ Check that the bases in the domains created for a structure are consistent with the bases in the strand
            they are part of.

            The combination of the bases in a list of domains for a strand should equal the number of base and follow the 
            order of bases in that strand. In addition each domain should only contain bases for a single helix. 
        """
        self._logger.setLevel(logging.INFO)
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("============================== check domains ============================== " )
        num_failures = 0
        for strand in self.strands:
            self._logger.debug("-------------------- strand %d -------------------- " % strand.id)
            domain_list = strand.get_domains()
            self._logger.debug("Number of domains: %d" % len(domain_list))

            # Get the combined domain bases and check that each domain is in only one helix.
            domain_base_ids = []
            for domain in domain_list:
                helix = domain.base_list[0].h
                helix_list = [ helix ]
                self._logger.debug("Domain %d: number of bases: %d" % (domain.id, len(domain.base_list)))
                for base in domain.base_list:
                    if base.h != helix:
                        helix = base.h
                        helix_list.append(helix)
                    domain_base_ids.append(base.id)
                if len(helix_list) != 1:
                    self._logger.error("The domain %d references more than one helix: %s" % str(helix_list))
            #_for domain in domain_list

            if len(strand.tour) != len(domain_base_ids):
                self._logger.error("The number of domain bases %d does not equal the number of strand bases %d." %
                    (len(domain_base_ids), len(strand.tour)))
                self._logger.error("Strand bases: %s" % (str(strand.tour)))
                self._logger.error("Domain bases: %s" % (str(domain_base_ids)))
                num_failures += 1
                continue 

            # Check that the combined domain bases are equal to the bases in the strand. 
            match_failed = False
            for sbid,dbid in zip(strand.tour,domain_base_ids):
                if sbid != dbid: 
                    match_failed = True
                    self._logger.error("The domain base %d does not match the strand base %d." % (dbid, sbid))
                    num_failures += 1
                    break
            #__for sbid,dbid in zip(strand.tour,domain_base_ids)

            # For a failed match check if domain bases are out of order.
            if match_failed:
                id_set = set()
                for id in strand.tour:
                    id_set.add(id)
                num_match_failed = 0
                for id in domain_base_ids:
                    if id not in id_set:
                        num_match_failed += 1
                if num_match_failed == 0:
                    self._logger.error("Check failed: Domain bases match strand bases but are not in the same order.")
            else:
                self._logger.debug("Check passed: domain bases match strand bases.")
            #__if match_failed
        #__for strand in self.strands

        if num_failures == 0:
            self._logger.info("Domain consistency check: all domains passed.")
        else:
            self._logger.error("Domain consistency check: %d domains failed." % num_failures)

    def _add_domain(self, id, strand, base_list, msg=""):
        """ Create a DnaDomain object from a list of bases. 

            Arguments:
                id (int): The domain ID; starts from 0.
                strand (DnaStrand): The strand the domain is in.
                base_list (list[DnaBase]): The list of bases in the domain.
                msg (string): An optional string used for debugging. 

            The new domain object is added to self.domain_list and to the list of domains for the helix and 
            strand it is contained in.
        """ 
        start_pos = base_list[0].p
        end_pos = base_list[-1].p
        self._logger.debug("++++ add domain %d   %s" % (id,msg))
        self._logger.debug("     vh:  %4d  -  %4d" % (base_list[0].h, base_list[-1].h))
        self._logger.debug("     pos: %4d  -  %4d" % (base_list[0].p, base_list[-1].p))
        base = base_list[0]
        helix = self.structure_helices_map[base.h]
        domain = nd.Domain(id, helix, strand, base_list)
        helix.domain_list.append(domain)
        strand.domain_list.append(domain)
        for base in base_list:
            base.domain = domain.id
        self.domain_list.append(domain)

    def _compute_strand_helix_references(self):
        """ Set the virtual helices referenced by each strand. """
        for strand in self.strands:
            for id in strand.tour:
                base = self.base_connectivity[id-1]
                helix = self.structure_helices_map[base.h]
                strand.add_helix(helix)
        #__for strand in self.strands__

    def _set_helix_connectivity(self):
        """ For each helix set the list of helices it is connected to. """ 
        self._logger.debug("[DnaModel::==================== set_vhelix_connectivity==================== ] ")
        for helix1 in self.structure_helices:
            self._logger.debug(" ----- vhelix num %d -----" % helix1.lattice_num)
            helix_connectivity = []
            row = helix1.lattice_row
            col = helix1.lattice_col 
            for helix2 in self.structure_helices:
                if helix1 == helix2: 
                    continue 
                if (abs(row-helix2.lattice_row) + abs(col-helix2.lattice_col) < 2) and \
                    self.lattice.get_neighbor_index(row, col, helix2.lattice_row, helix2.lattice_col) != -1:
                    connection = DnaHelixConnection(helix1,helix2)
                    helix_connectivity.append(connection)
                    self._logger.debug("connected to %d " % helix2.lattice_num)
            helix1.helix_connectivity = helix_connectivity

    def _compute_helix_design_crossovers(self):
        """ Compute the design cross-overs for all helices.
        """
        for helix in self.structure_helices:
            helix.compute_design_crossovers(self)

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
        self.helix_connectivity = []
        self.start_pos = -1
        self.start_staple = -1
        self.start_scaffold = -1
        self.possible_staple_crossovers = []
        self.possible_scaffold_crossovers = []

        # note [davep] are we using these?
        self.triad = None
        self.id_nt = None
        self.helix_topology = None

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

    def compute_design_crossovers(self,dna_structure):
        logger = self._setup_logging('DnaStructureHelix'+str(self.id))
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
                            crossover = DnaCrossover(self,connection,base,strand)
                            connection.crossovers.append(crossover)
                    #__if (down != -1)

                    if (up != -1):
                        up_base = base_connectivity[up-1]
                        if (up_base.h != base.h) and (up_base.h == num):
                            logger.debug("base:%4d  p:%4d  h:%4d" % (base.id, base.p, base.h))
                            logger.debug("  xu:%4d  p:%4d  h:%4d" % (up_base.id, up_base.p, up_base.h))
                            strand = dna_structure.get_strand(base.strand)
                            crossover = DnaCrossover(self,connection,base,strand)
                            connection.crossovers.append(crossover)
                    #__if (up != -1)

                #__for base in self.base_list
            #__for base_list in [self.staple_base_list,self.scaffold_base_list]
            logger.debug(">>> added %d crossovers " % len(connection.crossovers))
        #__for connection in self.helix_connectivity:

class DnaHelixConnection(object):
    """ This class stores information for connected helices.
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
        #print("[DnaHelixConnection] --------- compute direction ---------") 
        # Get the first helix axis and a point on that axis.
        helix1 = self.from_helix
        #print(">>> helix1: num: %d  row: %d  col: %d" % (helix1.lattice_num, helix1.lattice_row, helix1.lattice_col)) 
        start_pos = next((i for i in xrange(0,len(helix1.staple_base_list)) if helix1.staple_base_list[i] != None),-1)
        helix1_base = helix1.staple_base_list[start_pos]
        pt1 = helix1.helix_axis_nodes[helix1_base.p]
        axis1 = [helix1.end_frames[0,2,0], helix1.end_frames[1,2,0], helix1.end_frames[2,2,0]]

        # Get the second (adjacent) helix axis and a point on that axis.
        helix2 = self.to_helix
        #print(">>> helix2: num: %d  row: %d  col: %d" % (helix2.lattice_num, helix2.lattice_row, helix2.lattice_col))
        start_pos = next((i for i in xrange(0,len(helix2.staple_base_list)) if helix2.staple_base_list[i] != None),-1)
        helix2_base = helix2.staple_base_list[start_pos]
        pt2 = helix2.helix_axis_nodes[helix2_base.p]
        axis2 = [helix2.end_frames[0,2,0], helix2.end_frames[1,2,0], helix2.end_frames[2,2,0]]
        axis2_length = np.linalg.norm(axis2)

        # Compute the unit vector in the direction of the adjacent helix.
        vec = pt1 - pt2
        d = np.dot(axis2,vec) / axis2_length
        a2pt = pt2 + np.dot(axis2,d)
        self.direction = a2pt - pt1
        self.direction = self.direction / np.linalg.norm(self.direction)
        #print(">>> direction: %g %g %g" % (self.direction[0], self.direction[1], self.direction[2]))

class DnaCrossover(object):
    """ This class stores information for a cross-over between two helices.

        Attributes:
            helix (DnaStructureHelix): 
    """
    def __init__(self, helix, helix_connection, crossover_base, strand):
        self.helix = helix
        self.helix_connection = helix_connection
        self.crossover_base = crossover_base
        self.strand = strand




