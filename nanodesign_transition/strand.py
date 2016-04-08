#!/usr/bin/env python
"""
This module is used to store information for a DNA strand.

A DNA strand is a continuous chain of nucleotides. It can be either a scaffold of a staple.
"""
import sys
import os
from sets import Set
import numpy as np

class DnaStrand(object):
    def __init__(self, id):
        self.id = id
        self.is_scaffold = False
        self.is_circular = False
        self.tour = []
        self.is_main = []
        self.seq = []
        self.insert_seq = []
        self.rotations = []
        self.translations = []
        self.color = [1.0,1.0,1.0]
        self.helix_list = dict()
        self.base_id_list = dict()
        self.dna_structure = None
        self.base_coords = None

    def add_helix(self, helix): 
        id = helix.lattice_num
        if (id not in self.helix_list):
            #print("[DnaStrand] ---------- strand %d  add helix %d ----------" % (self.id, id))
            self.helix_list[id] = helix

    def get_domains_info(self):
        #print("")
        #print("[DnaStrand] ---------- get domains info ----------")
        #print("[DnaStrand] is_scaffold " + str(self.is_scaffold))
        domains_info = []
        domain_list = self.get_domains()
        for domain in domain_list:
           did = domain.id
           domains_info.append(did) 
           #print("[DnaStrand] domain id %d  info %s " % (id, domain.get_) 
        return domains_info 

    def get_domains(self):
        #print("[DnaStrand] ---------- get domains ----------")
        #print("[DnaStrand] self.helix_list size %d: " % len(self.helix_list)) 
        domain_list = dict()
        for id, helix in self.helix_list.iteritems():
           helix_domains = helix.domain_list
           #print("[DnaStrand] helix id %d   num domains %d " % (id,len(helix_domains))) 
           for domain in helix_domains:
               id = domain.id
               strand_id = domain.strand.id
               #print("[DnaStrand] domain id %d " % id) 
               #if ((id not in domain_list) and (strand_id == self.id)):
               if ((strand_id == self.id) and (id not in domain_list)): 
                   domain_list[id] = domain
                   #print("[DnaStrand] add domain id %d " % id) 
           #__for domain in domain_list
        #__for id, helix in self.helix_list.iteritems

        # Create a domain list sorted by strand position.
        domain_base_pos = {}
        for id in domain_list:
            domain = domain_list[id]
            base = domain.base_list[0]
            bindex = self.get_base_index(base)
            domain_base_pos[bindex] = domain
        #print("######## sorted domains ##########")
        domain_sorted_list = []
        for pos in sorted(domain_base_pos):
            domain = domain_base_pos[pos]
            domain_sorted_list.append(domain_base_pos[pos])
            #print(">>> domain id: %d" % (domain.id))
        return domain_sorted_list

    def get_base_coords(self):
        """ Get the coordinates of bases along the dna helix axis. """
        if (not self.base_coords):
            num_bases = len(self.tour)
            self.base_coords = np.zeros((num_bases,3), dtype=float)
            for i in xrange(0,num_bases):
                id = self.tour[i]
                base = self.dna_structure.base_connectivity[id-1]
                helix_num = base.h
                helix_pos = base.p
                helix = self.helix_list[helix_num]
                nodes = helix.helix_axis_nodes
                self.base_coords[i] = nodes[base.p]
        #__if (not self.base_coords)
        return self.base_coords

    def get_base_index(self, base):
        """ Get the index into the strand for the given base. """
        if (not self.base_id_list):
            num_bases = len(self.tour)
            for i in xrange(0,num_bases):
                id = self.tour[i]
                self.base_id_list[id] = i
        return self.base_id_list[base.id]


