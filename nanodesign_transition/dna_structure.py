#!/usr/bin/env python
"""
This module defines the classes used to define the connectivity and geometry of a DNA structure. 

A DNA structure consists of a number of scaffold and staple strands (DNA origami), or oligo strands alone, bound together 
to form a designed geometric shape.
"""
from collections import OrderedDict
import json
import logging
import numpy as np
from parameters import DnaParameters
from converters.cadnano.common import CadnanoLatticeType
from dna_structure_helix import DnaStructureHelix,DnaHelixConnection

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

class DnaStructure(object):
    """ This class stores the base connectivity and geometry for a DNA model. 

        Attributes:
            base_connectivity (List[DnaBase]): The list of DnaBase objects for the structure.
            domain_list (List[Domain]): The list of Domain objects for the structure.
            id_nt (NumPy Nx2 ndarray[int]): The base IDs for scaffold bases and their paired staple base.
            lattice_type (CadnanoLatticeType): The lattice type the geometry of this structure is derived from. 
            lattice (Lattice): The Lattice object used for calculating lattice-dependent data (e.g. neighboring lattice
                locations).  
            parameters (DnaParameters): Stores information for DNA parameters (e.g. helix radius).
            strands (List[DnaStrand]): The list a DnaStrand objects. 
            strands_map (Dict[DnaStrand]): The dictionary that maps strand IDs to DnaStrand objects.
    """ 

    def __init__(self, name, base_connectivity, helices, dna_parameters):
        """ Initialize a DnaStructure object. 

            Arguments:
                name (string): The name of the structure.
                base_connectivity (List[DnaBase]): The list of DNA bases for the structure. 
                helices (List[DnaStructureHelix]): The list of helices for the structure. 
                dna_parameters (DnaParameters): The DNA parameters to use when creating the 3D geometry for the design.
        """ 
        self.name = name
        self.lattice_type = CadnanoLatticeType.none
        self.lattice = None
        self.dna_parameters = dna_parameters
        self.base_connectivity = base_connectivity
        self.structure_helices_map = dict()
        self.structure_helices_coord_map = dict()
        self.strands = None
        self.strands_map = dict()
        self.domain_list = []
        self.connector_points = []
        self._logger = self._setup_logging()
        self._add_structure_helices(helices)

    def _setup_logging(self):
        """ Set up logging."""
        logger = logging.getLogger(__name__ + ":" + str(self.name))
        logger.setLevel(logging.INFO)
        # Create console handler and set format.
        if not len(logger.handlers):
            console_handler = logging.StreamHandler()
            formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
            console_handler.setFormatter(formatter)
            logger.addHandler(console_handler)
        return logger

    def _add_structure_helices(self, structure_helices):
        """ Add a list of structural helices. 

            The structural helices are stored in two dictionaries, one used to look up a helix 
            using a helix ID and the other used to look up a helix using a lattice (row,col) 
            coordinate.
            # TODO (DaveP) We need to remove references to lattice-based geometry.

            Arguments:
               structure_helices (list[DnaStructureHelix]): The list of structural helices for the structure.
        """
        for helix in structure_helices:
            self.structure_helices_map[helix.id] = helix
            self.structure_helices_coord_map[(helix.lattice_row,helix.lattice_col)] = helix

    def set_lattice_type(self, lattice_type):
        """ Set the lattice type the geometry of this structure is derived from. 

            Arguments:
                lattice_type(CadnanoLatticeType): The type of lattice (e.g. square, honeycomb). 

            A Lattice object is created from the given lattice type. It is used to calculate 
            lattice-dependent data such as crossover locations and the neighboring lattice 
            coordinates for a given coordinate.
        """
        self.lattice_type = lattice_type
        self.lattice = Lattice.create_lattice(lattice_type, self.dna_parameters.helix_distance)

    def compute_aux_data(self):
        """ Compute auxiliary data. 

            Compute data derived from the base connectivty: domains, strand/helix and helix/helix connectivty 
            relationships, and crossovers. This data is needed for visualization, calculating melting
            temperature and other applications.
        """
        for strand in self.strands:
            strand.dna_structure = self
        self.set_strand_helix_references()
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

    def set_strand_helix_references(self):
        """ Set the helices referenced by each strand. """
        for strand in self.strands:
            for base in strand.tour:
                helix = self.structure_helices_map[base.h]
                strand.add_helix(helix)
        #__for strand in self.strands__

    def apply_helix_xforms(self, helix_group_xforms):
        """ Apply helix group transformations.

            Arguments:
                helix_group_xforms (List[HelixGroupXform]): The list of helix group transforms.

            The geometry for the list of helices for each group are rotated and translated together 
            by the given transformation. 
        """
        self._logger.setLevel(logging.INFO)
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("===================== apply helix group xforms =====================")
        for helix_group in helix_group_xforms:
            # Get the center of the helix group.
            group_center = np.array([0.0,0.0,0.0], dtype=float)
            for helix in helix_group.helices:
                helix_center = helix.get_center()
                group_center += helix_center 
            #__for helix_id in helix_group.helices
            group_center /= len(helix_group.helices)
            helix_group.transformation.set_center(group_center)
 
            # Transform helices geometry.
            for helix in helix_group.helices:
                helix.apply_xform(helix_group.transformation)
            #__for helix_id in helix_group_xforms.helices:
        #__for helix_group in helix_group_xforms

    def xform_from_connectors(self, connector_strands, helix_ids, xform):
        """ Create a transformation to rotate and translate a group of helices so that
            distance crossovers (connections) align.

            Arguments:
                connector_strands (List[DnaStrand]): The list of strands that contain connections. 
                helix_ids (List[int]): The list of helix IDs to transform.
                xform (Xform): The transformation to set the rotation and translation used to 
                    align the crossovers.

            This is an attempt to automatically generate a 3D structure from distanct scaffold
            crossovers used to assemble groups of helices. The algorithm builds two lists of 
            points for the scaffold connections and attempts to align them. This can typically produce 
            a transformation that only translates the two sections together but does not oriented them,  
            the relative rotation of two components created by excluded volume interactions. 
        """
        self._logger.setLevel(logging.INFO)
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("================ xform_from_connectors ==================") 
        # Create a set of helix IDs to include in the transformation.
        helix_set = set()
        for id in helix_ids:
            helix_set.add(id)

        # Create the points lists for the connections. 
        max_dist = 2.0*self.dna_parameters.helix_distance
        points1 = []
        points1_map = set()
        points2 = []
        points2_map = set()
        for strand in connector_strands:
            tour = strand.tour 
            for i in xrange(0,len(tour)-1):
                base1 = tour[i]
                pt1 = base1.coordinates
                base2 = tour[i+1]
                pt2 = base2.coordinates
                dist = np.linalg.norm(pt1 - pt2)
                if dist > max_dist:
                    self._logger.debug(" ")
                    self._logger.debug("Base 1  h %d  p %d  (%g %g %g) " % (base1.h, base1.p, 
                        base1.coordinates[0], base1.coordinates[1], base1.coordinates[2] ))
                    self._logger.debug("   down h %d  p %d  " % (base1.down.h, base1.down.p))
                    self._logger.debug("   up   h %d  p %d  " % (base1.up.h, base1.up.p))

                    self._logger.debug("Base 2  h %d  p %d  (%g %g %g) " % (base2.h, base2.p, 
                        base2.coordinates[0], base2.coordinates[1], base2.coordinates[2] ))

                    if (base1.h,base1.p) not in points1_map:
                        points1_map.add((base1.h,base1.p))
                        if base1.h in helix_set:
                            points1.append(pt1)
                            self._logger.debug("Add Base 1 to points1")
                        else:
                            points2.append(pt1)
                            self._logger.debug("Add Base 1 to points2")
                    #__if (base1.h,base1.p) not in points1_map

                    if (base2.h,base2.p) not in points2_map:
                        points2_map.add((base2.h,base2.p))

                        if base2.h in helix_set:
                            points1.append(pt2)
                            self._logger.debug("Add Base 2 to points1")
                        else:
                            points2.append(pt2)
                            self._logger.debug("Add Base 2 to points2")
                    #__if (base2.h,base2.p) not in points2_map:
            #__for i in xrange(0,len(tour)-1)
        #__for strand in connector_strands
        self._logger.debug("Number of connector points %d" % len(points1))

        # Set the transformation that fits points1 to points2.
        xform.rms_fit(points1, points2)

        # Set the connector points for later visualization.
        self.connector_points = [points1, points2]

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

        # Set flag for merging domains.
        merge_domains = True
        merge_domains = False

        # Iterate over the scaffold and staple strands of a structure. 
        for strand in self.strands:
            self._logger.debug("")
            if ( strand.is_scaffold):
                self._logger.debug("==================== scaffold strand %d ====================" % strand.id)
            else:
                self._logger.debug("==================== staple strand %d ====================" % strand.id)

            start_base = strand.tour[0]
            end_base = strand.tour[-1]
            self._logger.debug("Strand number of bases: %3d" % len(strand.tour))
            self._logger.debug("Strand start: h: %3d  p: %3d" % (start_base.h, start_base.p))
            self._logger.debug("Strand end: h: %3d  p: %3d" % (end_base.h, end_base.p))

            # Initialize the domain base list.
            base = strand.tour[0]
            curr_across_sign = 0 if base.across else -1
            domain_bases = [ base ]

            # Traverse the bases in a strand and create domains.
            for i in xrange(1,len(strand.tour)):
                base = strand.tour[i]
                across_sign = 0 if base.across else -1
                add_curr_base = True
                #self._logger.debug("Base h %3d  p %3d " % (base.h, base.p))
                #if (base.up): 
                #    self._logger.debug("    Base up  h %3d  p %3d " % (base.up.h, base.up.p))

                # If no domain bases then just continue after checking for sign change.
                if len(domain_bases) == 0:
                    domain_bases.append(base)
                    if curr_across_sign != across_sign:
                        curr_across_sign = across_sign
                    continue

                # Check for a single->double or double->single strand transition.
                if curr_across_sign != across_sign:
                    domain_id = self._add_domain(domain_id, strand, domain_bases, merge_domains, "sign change")
                    domain_bases = []
                    curr_across_sign = across_sign

                # Check for a crossover between helices for this base.
                elif self._check_base_crossover(base):
                    last_base = domain_bases[-1]
                    # Make sure the current base is in the same helix.
                    if base.h == last_base.h:
                        domain_bases.append(base)
                        add_curr_base = False
                    domain_id = self._add_domain(domain_id, strand, domain_bases, merge_domains, "base crossover")
                    domain_bases = []

                # Check the base paired to this base for: crossover or termination.
                elif base.across != None:
                    abase = base.across
                    if self._check_base_crossover(abase):
                        domain_bases.append(base)
                        add_curr_base = False
                        domain_id = self._add_domain(domain_id, strand, domain_bases, merge_domains, "abase crossover")
                        domain_bases = []
                    # If a strand terminates make sure the current base is in the same strand. 
                    elif (abase.down == None) or (abase.up == None):
                        last_base = domain_bases[-1]
                        if last_base.across != None:
                            last_abase = last_base.across
                            if abase.strand == last_abase.strand:
                                domain_bases.append(base)
                                add_curr_base = False
                        else:
                            domain_bases.append(base)
                            add_curr_base = False
                        domain_id = self._add_domain(domain_id, strand, domain_bases, merge_domains, "abase start/end")
                        domain_bases = []
                #__if curr_across_sign != across_sign

                # Add the current base to the current list of domain bases.
                if add_curr_base:
                    domain_bases.append(base)
            #__for i in xrange(1,len(strand.tour))

            # Add a domain for any remaining bases.
            if len(domain_bases) != 0:
                domain_id = self._add_domain(domain_id, strand, domain_bases, merge_domains, "remaining")
        #__for strand in self.strands__

        self._logger.info("Number of domains computed: %d " % len(self.domain_list))

        # Check if the computed domains are consistent with the strands they were computed from.
        self.check_domains()

        # Set the strand and domain each domain is connected to.
        for domain in self.domain_list:
            across_base = None
            for base in domain.base_list:
                if (base.across != None):
                    across_base = base.across
                    break
            #__for base in domain.base_list
            conn_dom = -1
            conn_strand = -1
            if (across_base != None):
                conn_dom = across_base.domain
                conn_strand = across_base.strand
            domain.connected_strand = conn_strand
            domain.connected_domain = conn_dom
        #__for domain in domain_list

    def _check_base_crossover(self, base):
        """ Check if there is a crossover to a different helix at the given base.

            Arguments:
                base (DnaBase): The base to check for a crossover.

            Returns True if there is a crossover at the base.
        """
        if base.down == None:
            return False
        if base.down.h != base.h:
            return True

        if base.up == None:
            return False
        if base.up.h != base.h:
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
            self._logger.debug("Number of bases %d " % len(strand.tour))
            strand_bases = ""
            for base in strand.tour:
                strand_bases += " " + str(base.id)
            self._logger.debug("Bases: %s " % strand_bases) 
            domain_list = strand.domain_list
            self._logger.debug("Number of domains: %d" % len(domain_list))

            # Get the combined domain bases and check that each domain is in only one helix.
            domain_base_ids = []
            for domain in domain_list:
                helix = domain.base_list[0].h
                helix_list = [ helix ]
                self._logger.debug("Domain %d: number of bases: %d" % (domain.id, len(domain.base_list)))
                for base in domain.base_list:
                    self._logger.debug("       base id %d  h %d  p %d" % (base.id, base.h, base.p))
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
                self._logger.error("Strand bases: %s" % (str(strand_bases)))
                self._logger.error("Domain bases: %s" % (str(domain_base_ids)))
                num_failures += 1
                continue 

            # Check that the combined domain bases are equal to the bases in the strand. 
            match_failed = False
            if (not strand.is_circular) :
                for sbase,dbid in zip(strand.tour,domain_base_ids):
                    if sbase.id != dbid: 
                        match_failed = True
                        self._logger.error("The domain base %d does not match the strand base %d." % (dbid, sbid))
                        num_failures += 1
                        break
                #__for sbase,dbid in zip(strand.tour,domain_base_ids)
            #__if (not strand.is_circular)

            # For a failed match check if domain bases are out of order.
            if match_failed:
                id_set = set()
                for base in strand.tour:
                    id_set.add(base.id)
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



    def _add_domain(self, id, strand, base_list, merge_domains, msg=""):
        """ Create a DnaDomain object from a list of bases. 

            The new domain object is added to self.domain_list and to the list of domains for the helix and 
            strand it is contained in.

            If a strand is circular then the first domain will be created from the bases at the start of the strand. 
            The remaining bases at the end of the strand will be merged into the initial domain. 

            The merging algorithm works as follows:

                first_dom = the first domain created for the strand.
                first_sbase = the first base in the first domains base list.
                base_list = list of bases for a new domain 
                ebase = the last base in base_list

                      helix 1     .------------se-------------------------.
                                  |                                       |
                      helix 2     .---------------------------------------.


                                          first_sbase
                                               |
                      1st domain [.------------s]                            

                      new domain                [e-------------------------.]
                                                 |
                                                 ebase

                if first_base and ebase are in the same helix and they are a sigle base apart then add the new domains
                base list to the first domains:

                      1st domain  [.------------se-------------------------.]


            TODO (DaveP) We need to check for the edge case where a strand crosses over to another helix
                         a single base away from the strand start. 

                                    .----------.
                                    |          |
                                    .---..-----.               s = strand start
                                        ||                     . = crossover 
                                  .-----.s---------.
                                  |                |
                                  .----------------.
             
             

            Arguments:
                id (int): The domain ID; starts from 0.
                strand (DnaStrand): The strand the domain is in.
                base_list (list[DnaBase]): The list of bases in the domain.
                merge_domains (bool): If True then for circular strands merge the bases from the start of the strand with
                    those from the end. 
                msg (string): An optional string used for debugging. 

            Returns:
                id (int): The next domain ID. If a domain is added then id is incremented by one and returned.

        """ 
        # Check if the bases should be merged into the strand's first domain.
        domain_was_merged = False
        if merge_domains and strand.is_circular and len(strand.domain_list):
            first_dom = strand.domain_list[0]
            first_sbase = first_dom.base_list[0]
            first_ebase = first_dom.base_list[-1]
            sbase = base_list[0]
            ebase = base_list[-1]

            if (first_sbase.h == ebase.h) and  (abs(first_sbase.p-ebase.p) == 1):
                #self._logger.setLevel(logging.DEBUG)
                self._logger.debug("---------- strand %d: merging domains ---------- " % strand.id)
                self._logger.debug(">>> first domain id %d:  h %d  sb %d  eb %d" % (first_dom.id, first_sbase.h, first_sbase.p, first_ebase.p))
                self._logger.debug(">>> base list:           sbase %d  ebase %d" % (sbase.p, ebase.p))
                self._logger.setLevel(logging.INFO)
                domain_was_merged = True 
                merged_base_list = []
                for base in base_list:
                    base.domain = first_dom.id
                    merged_base_list.append(base)
                for base in first_dom.base_list:
                    merged_base_list.append(base)
                first_dom.base_list = merged_base_list
        #__if strand.is_circular and len(strand.domain_list):

        # If we haven't merged the bases then create a new domain.
        if not domain_was_merged:
            start_pos = base_list[0].p
            end_pos = base_list[-1].p
            self._logger.debug("++++ add domain %d   %s" % (id,msg))
            self._logger.debug("     vh:  %4d  -  %4d" % (base_list[0].h, base_list[-1].h))
            self._logger.debug("     pos: %4d  -  %4d" % (base_list[0].p, base_list[-1].p))
            base = base_list[0]
            helix = self.structure_helices_map[base.h]
            domain = nd.Domain(id, helix, strand, base_list)
            strand.domain_list.append(domain)
            for base in base_list:
                base.domain = domain.id
            self.domain_list.append(domain)
            id = id + 1
        #__if not domain_was_merged:
        return id

    def _compute_strand_helix_references(self):
        """ Set the virtual helices referenced by each strand. """
        for strand in self.strands:
            for base in strand.tour:
                helix = self.structure_helices_map[base.h]
                strand.add_helix(helix)
        #__for strand in self.strands__

    def _set_helix_connectivity(self):
        """ For each helix set the list of helices it is connected to. """ 
        self._logger.debug("[DnaModel::==================== set_vhelix_connectivity==================== ] ")
        for helix1 in self.structure_helices_map.itervalues():
            self._logger.debug(" ----- vhelix num %d -----" % helix1.lattice_num)
            helix_connectivity = []
            row = helix1.lattice_row
            col = helix1.lattice_col 
            for helix2 in self.structure_helices_map.itervalues():
                if helix1 == helix2: 
                    continue 
                if (abs(row-helix2.lattice_row) + abs(col-helix2.lattice_col) < 2) and \
                    self.lattice.get_neighbor_index(row, col, helix2.lattice_row, helix2.lattice_col) != None:
                    connection = DnaHelixConnection(helix1,helix2)
                    helix_connectivity.append(connection)
                    self._logger.debug("connected to %d " % helix2.lattice_num)
            helix1.helix_connectivity = helix_connectivity

    def _compute_helix_design_crossovers(self):
        """ Compute the design cross-overs for all helices.
        """
        for helix in self.structure_helices_map.itervalues():
            helix.compute_design_crossovers(self)

    def write(self, file_name, write_json_format):
        """ Write the structure information to a file. 
            Structure information is written to files in JSON and plain text formats.
        """

        # Compute auxillary data to calculate domains.
        self.compute_aux_data()

        # Write structure information in JSON format.
        if write_json_format:
            self._logger.info("Writing DNA strcuture to file %s." % file_name)
            base_list = []
            for base in self.base_connectivity:
                base_info = OrderedDict()
                base_info['id'] = base.id
                base_info['helix'] = base.h
                base_info['pos'] = base.p
                base_info['up'] =  base.up.id if base.up else -1 
                base_info['down'] = base.down.id if base.down else -1
                base_info['across'] = base.across.id if base.across else -1
                base_info['sequence'] = base.seq
                base_info['strand'] = base.strand
                base_list.append(base_info)
            #__for base in self.base_connectivity

            strand_list = []
            for strand in self.strands:
                strand_info = OrderedDict()
                strand_info['id'] = strand.id
                strand_info['scaffold'] = strand.is_scaffold
                strand_info['bases'] = [base.id for base in strand.tour]
                strand_info['domain_ids'] = [domain.id for domain in strand.domain_list ]
                strand_list.append(strand_info)

            domain_list = []
            for domain in self.domain_list:
                domain_info = OrderedDict()
                domain_info['id'] = domain.id
                domain_info['bases'] = [base.id for base in domain.base_list]
                domain_list.append(domain_info)

            structure = OrderedDict()
            structure['num_bases'] = len(base_list)
            structure['num_strands'] = len(strand_list)
            structure['num_domains'] = len(domain_list)
            structure['bases'] = base_list
            structure['strands'] = strand_list
            structure['domains'] = domain_list 

            with open(file_name, 'w') as outfile:
                json.dump(structure, outfile, indent=4, separators=(',', ': '), sort_keys=False)

        # Write structure information in plain text format.
        file_name = file_name.replace("json", "txt")
        self._logger.info("Writing DNA structure in text format to file %s." % file_name)
        with open(file_name, 'w') as outfile:
            outfile.write("# number of bases %d\n" % len(self.base_connectivity))
            outfile.write("# number of strands %d\n" % len(self.strands))
            outfile.write("# number of domains %d\n" % len(self.domain_list))
            outfile.write("# bases: id   helix  pos   up   down  across  seq   strand   scaf\n")
            for base in self.base_connectivity:
                up =  base.up.id if base.up else -1 
                down = base.down.id if base.down else -1
                across = base.across.id if base.across else -1
                outfile.write("%4d %5d %5d %5d %5d %5d  %5s  %5d  %5d\n" % \
                    (base.id, base.h, base.p, up, down, across, base.seq, base.strand, base.is_scaf))
            outfile.write("# strands: id  scaf  numBases:[baseIDs]  numDomains:[domainIDs]\n")
            for strand in self.strands:
                outfile.write("strand %4d %2d \n" % (strand.id, strand.is_scaffold))
                outfile.write("%d:%s\n" % (len(strand.tour), str([base.id for base in strand.tour])))
                outfile.write("%d:%s\n" % (len(strand.domain_list), [domain.id for domain in strand.domain_list ]))
            outfile.write("# domains: id  numBases:[baseIDs]\n")
            for domain in self.domain_list:
                outfile.write("domain %4d   %d:%s\n" % \
                    (domain.id, len(domain.base_list), [base.id for base in domain.base_list]))
        #__with open(file_name, 'w') as outfile

    def write_topology(self, file_name, write_json_format):
        """ Write the base information with base connectivity to a file. 
            Base information is written to files in JSON and plain text formats.
        """
        # Write base information in JSON format.
        if write_json_format:
            self._logger.info("Writing DNA base connectivity in JSON format to file %s." % file_name)
            base_list = []
            for base in self.base_connectivity:
                base_info = OrderedDict()
                base_info['id'] = base.id 
                base_info['helix'] = base.h 
                base_info['pos'] = base.p 
                base_info['up'] =  base.up.id if base.up else -1 
                base_info['down'] = base.down.id if base.down else -1
                base_info['across'] = base.across.id if base.across else -1
                base_info['sequence'] = base.seq 
                base_info['strand'] = base.strand 
                base_list.append(base_info)
            #__for base in self.base_connectivity

            topology = { 'bases' : base_list } 
            with open(file_name, 'w') as outfile:
                json.dump(topology, outfile, indent=4, separators=(',', ': '), sort_keys=False)

        # Write base information in plain text format.
        file_name = file_name.replace("json", "txt")
        self._logger.info("Writing DNA base connectivity in plain text format to file %s." % file_name)
        with open(file_name, 'w') as outfile:
            outfile.write("# id   helix  pos   up   down  across  seq   strand   scaf\n")
            for base in self.base_connectivity:
                up =  base.up.id if base.up else -1 
                down = base.down.id if base.down else -1
                across = base.across.id if base.across else -1
                outfile.write("%4d %5d %5d %5d %5d %5d  %5s  %5d  %5d\n" % \
                    (base.id, base.h, base.p, up, down, across, base.seq, base.strand, base.is_scaf))
        #__with open(file_name, 'w') as outfile

#__class DnaStructure(object):

