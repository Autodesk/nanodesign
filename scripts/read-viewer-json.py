#!/usr/bin/python
# -*- coding: utf-8 -*-
#================================================================================#
#                                 read a vis json file                           #
#================================================================================#
import json
import os
import sys

class VisJsonFields:
    LATTICE_TYPE = "lattice_type"
    STRANDS = "strands"
    BASES   = "bases"
    COLOR   = "color"
    DOMAINS = "domains"
    STRAND_ID = "strand_id"
    NUMBER_OF_BASES = "number_of_bases"
    VIRTUAL_HELICES = "virtual_helices"
    IS_SCAFFLOD = "is_scaffold"

#------------------------------------------------------------#
#                       main                                 #
#------------------------------------------------------------#
def main():

    if (len(sys.argv) == 2):
        file_name = sys.argv[1]
    else:
        file_name = "fourhelix_viewer.json"

    with open(file_name) as file:
        json_data = json.load(file)

    lattice_type = json_data[VisJsonFields.LATTICE_TYPE]
    print(">>> lattice type %s" % lattice_type) 

    print("========================= virtual helices =========================") 
    vhelix_list = json_data[VisJsonFields.VIRTUAL_HELICES]
    print(">>> number of virtual helices %d" % len(vhelix_list)) 
    for vhelix in vhelix_list:
        id = vhelix["id"]
        cadnano_info = vhelix["cadnano_info"]
        num = cadnano_info["num"]
        print("")
        print("--------------- vhelix %s --------------" % id)
        domain_list = vhelix[VisJsonFields.DOMAINS]
        print(">>> num %d" % num)
        print(">>> number of domains %d" % len(domain_list)) 
        print(">>> domains %s" % str(domain_list)) 
        strand_radius = vhelix["strand_radius"]
        print(">>> strand_radius %f" % strand_radius) 

    print("")
    print("========================= strands =========================") 
    strand_list = json_data[VisJsonFields.STRANDS]
    print(">>> number of strands %d" % len(strand_list)) 
    for strand in strand_list:
        id = strand["id"]
        is_scaffold = strand[VisJsonFields.IS_SCAFFLOD]
        bases = strand[VisJsonFields.BASES]
        vhelix_list = strand[VisJsonFields.VIRTUAL_HELICES]
        domain_list = strand[VisJsonFields.DOMAINS]
        color = strand[VisJsonFields.COLOR]
        print("")
        print(">>> strand %s: scaffold %s nbases %d nvhelix %d ndoms %d doms %s" % (id, is_scaffold, len(bases), len(vhelix_list),
           len(domain_list), str(domain_list) ))
        print("              color %s" % str(color)) 
        print("            : bases : ") 
        for base in bases:
            id = base['id']
            coord = base['coordinates']
            seq = base['sequence']
            print("            : id %d  coord %s  seq %s" % (id, str(coord), seq))

    print("")
    print("========================= domains =========================") 
    domains_list = json_data[VisJsonFields.DOMAINS]
    print(">>> number of domains %d" % len(domains_list)) 
    for domain in domains_list:
        id = domain["id"]
        num_bases = domain["number_of_bases"]
        base_list = domain["bases"]
        strand_id = domain["strand_id"]
        connected_strand = domain['connected_strand']
        connected_domain = domain['connected_domain']
        start_base_index = domain['start_base_index']
        end_base_index = domain['end_base_index']
        orientation = domain['orientation']
        start_pos = domain['start_position']
        end_pos  = domain['end_position']
        color = domain[VisJsonFields.COLOR]

        #print("")
        print(">>> domain %s: strand %d  nbases %d  bases %s cstrand %d cdomain %d " % 
           (id, strand_id, num_bases, str(base_list), connected_strand, connected_domain))
        print("               start_pos (%f %f %f) end_pos (%f %f %f) " % (start_pos[0], start_pos[1], start_pos[2], 
            end_pos[0], end_pos[1], end_pos[2]))
        print("               orientation (%f %f %f) " % (orientation[0], orientation[1], orientation[2]))
        print("               start_base_index %d  end_base_index %d " % (start_base_index, end_base_index))
        print("               color %s " % (str(color)))

if __name__=="__main__":
    main()
    

