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
        self._compute_strand_helix_references()
        self._compute_domains()

    def get_domains(self):
        if (not self.domain_list): 
            self._compute_domains()
        return self.domain_list

    def _set_strand_helix_references(self):
        """ Set the helices referenced by each strand. """
        for strand in self.strands:
            for id in strand.tour:
                base = self.base_connectivity[id-1]
                helix = self.structure_helices_map[base.h]
                strand.add_helix(helix)
        #__for strand in self.strands__

    def _compute_domains(self):
        """ Compute dna domains, contiguous sequences of nucleic acids separated by crossover junctions."""
        num_domains = 0
        self.domain_list = []
        debug = True
        debug = False

        for strand in self.strands:
            if debug: print("")
            if (strand.is_scaffold):
                #continue 
                scaffold = True
                if debug: print("---------- scaffold strand %d ----------" % strand.id)
            else:
                #continue 
                scaffold = False
                if debug: print("---------- staple strand %d ----------" % strand.id)
            color = strand.color
            id = strand.tour[0]
            base = self.base_connectivity[id-1]
            curr_helix_id = base.h
            helix = self.structure_helices_map[curr_helix_id]

            across_base = base.across
            curr_across_base_sign = np.sign(base.across)
            if ( base.across > 0):
                across_base = self.base_connectivity[base.across-1]
                curr_across_strand_id = across_base.strand
            else:
                curr_across_strand_id = -1

            if debug: 
                print("++++++ add initial domain %d, helix %d ++++++" % (num_domains,curr_helix_id))
                print(">>> base id %d  up %d  down %d  across %d  helix %d  pos %d  strand %d" % 
                    (id, base.up, base.down, base.across, base.h, base.p, base.strand))
            domain = nd.Domain(num_domains,helix)
            num_domains += 1
            domain.base_list.append(base)
            domain.color = color
            domain.strand = strand
            self.domain_list.append(domain)
            helix.domain_list.append(domain)

            for i in xrange(1,len(strand.tour)):
                id = strand.tour[i]
                base = self.base_connectivity[id-1]
                across_base_sign = np.sign(base.across)
                domain_added = False
                if debug: 
                    print(">>> base id %d  up %d  down %d  across %d  helix %d  pos %d  strand %d" % 
                          (id, base.up, base.down, base.across, base.h, base.p, base.strand))

                if (scaffold and (base.across > 0)):
                    across_base = self.base_connectivity[base.across-1]
                    if debug: 
                        print(">>> across base %d  helix %d  pos %d  strd %d" %
                          (across_base.id, across_base.h, across_base.p, across_base.strand))
                        print("")

                    if ((curr_across_strand_id != -1) and (curr_across_strand_id != across_base.strand )):
                        curr_helix_id = across_base.h
                        helix = self.structure_helices_map[curr_helix_id]
                        if debug: 
                            print("")
                            print("++++++ add new domain %d, helix %d ++++++" % (num_domains,curr_helix_id))
                            print("")
                        curr_across_strand_id = across_base.strand
                        domain_added = True
                        domain = nd.Domain(num_domains,helix)
                        num_domains += 1
                        domain.color = color
                        domain.strand = strand
                        self.domain_list.append(domain)
                        helix.domain_list.append(domain)
                        domain.base_list.append(base)
                        base.domain = domain.id
                    elif ( base.across > 0):
                        curr_across_strand_id = across_base.strand
                #__if (scaffold and (base.across > 0))

                if ((base.h != curr_helix_id) or (across_base_sign != curr_across_base_sign)): 
                    curr_helix_id = base.h
                    curr_across_base_sign = np.sign(base.across)
                    helix = self.structure_helices_map[curr_helix_id]

                    if debug: print("++++++ add new domain %d  helix %d ++++++" % (num_domains,base.h))

                    domain_added = True
                    domain = nd.Domain(num_domains,helix)
                    num_domains += 1
                    domain.color = color
                    domain.strand = strand
                    self.domain_list.append(domain)
                    helix.domain_list.append(domain)
                    domain.base_list.append(base)
                    base.domain = domain.id

                if (not domain_added):
                    domain.base_list.append(base)
                    base.domain = domain.id
            #__for i in xrange(1,len(strand.tour))__
        #__for strand in self.strands__

        # set the strand and domain connectivity for the domains
        if debug: print("")
        if debug: print(" Number of domains added %d " % (num_domains))

        for domain in self.domain_list:
            across = -1
            for base in domain.base_list:
                if (base.across != -1):
                    across = base.across
            #__for base in domain.base_list
            conn_dom = -1
            conn_strand = -1
            if (across != -1):
                across_base = self.base_connectivity[across-1]
                conn_dom = across_base.domain
                conn_strand = across_base.strand
            domain.connected_strand = conn_strand
            domain.connected_domain = conn_dom
            if debug:
                print("")
                print(">>> domain %d: strand %d helix %d  nbases %d  base across: %d  conn strand/dom %d %d" % 
                    (domain.id, domain.strand.id, domain.helix.id, len(domain.base_list), across, conn_strand,
                     conn_dom ))

        #sys.exit(0)
    #__def _compute_domains(self)

    #--------------------------------------------------
    # set the virtual helices referenced by each strand
    #--------------------------------------------------
    def _compute_strand_helix_references(self):
        for strand in self.strands:
            for id in strand.tour:
                base = self.base_connectivity[id-1]
                helix = self.structure_helices_map[base.h]
                strand.add_helix(helix)
        #__for strand in self.strands__


#__class DnaStructure(object):


class DnaStructureHelix(object):
    """ This class stores information for a DNA structure helix. 

        A structure helix is a region in a DNA structure that forms a cylindrical structural element.
        composed of two DNA helices.

        Attributes:
            id (int): Helix ID (1 - number of helices in a structure).
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

        # note [davep] are we using these?
        self.start_pos = -1
        self.triad = None
        self.id_nt = None
        self.helix_topology = None



