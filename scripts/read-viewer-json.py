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
    HELIX_CONNECTIVITY = "helix_connectivity"
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
        file_name = "./results/fourhelix_viewer.json"

    with open(file_name) as file:
        json_data = json.load(file)

    lattice_type = json_data[VisJsonFields.LATTICE_TYPE]
    print(">>> lattice type %s" % lattice_type) 

    print("========================= virtual helices =========================") 
    vhelix_list = json_data[VisJsonFields.VIRTUAL_HELICES]
    print(">>> number of virtual helices %d" % len(vhelix_list)) 
    domains_obj_list = json_data[VisJsonFields.DOMAINS]
    for vhelix in vhelix_list:
        id = vhelix["id"]
        cadnano_info = vhelix["cadnano_info"]
        num = cadnano_info["num"]
        print("")
        print("--------------- vhelix %s --------------" % id)
        domain_list = vhelix[VisJsonFields.DOMAINS]
        print(">>> cadnano vhelix num %d" % num)
        print(">>> number of domains %d" % len(domain_list)) 
        print(">>> domain IDs %s" % str(domain_list)) 
        print(">>> domains: ") 
        for domain_id in domain_list:
            for domain in domains_obj_list:
                id = domain["id"]
                if (domain_id == id):
                    num_bases = domain["number_of_bases"]
                    connected_domain = domain['connected_domain']
                    strand_id = domain["strand_id"]
                    print("id: %4d  nbases: %3d  strand: %3d  connected to domain: %4d" % (domain_id, num_bases, strand_id, connected_domain))
        helix_connectivity = vhelix[VisJsonFields.HELIX_CONNECTIVITY]
        num_possible_staple_crossovers = vhelix['num_possible_staple_crossovers']
        num_possible_scaffold_crossovers = vhelix['num_possible_scaffold_crossovers']
        num_conn = 0
        for connection in helix_connectivity:
           if connection: 
               num_conn += 1
        print(">>> number of possible staple crossovers: %d" % num_possible_staple_crossovers) 
        print(">>> number of possible scaffold crossovers: %d" % num_possible_scaffold_crossovers) 
        print(">>> number of connected helices: %d" % num_conn) 
        print(">>> helix connections: ") 
        for connection in helix_connectivity:
           if not connection: 
               continue 
           helix_id = connection['helix_id']
           helix_num = connection['helix_num']
           angle = connection['angle']
           dir = connection['direction']
           crossovers = connection['crossovers']
           print("helix id: %d  num: %d  dir: (%g %g %g)  angle: %g" % (helix_id, helix_num, dir[0], dir[1], dir[2], angle))
           if not crossovers: 
               print("    crossovers: None") 
               continue 
           print("    number of crossovers: %d" % len(crossovers)) 
           #print("    crossovers data: %s" % str(crossovers)) 
           print("    crossovers:")
           n = 1
           for crossover in crossovers:
               vhelix_base_index = crossover['vhelix_base_index']
               strand1_ID = crossover['first_strand_ID']
               strand1_bindex = crossover['first_strand_base_index']
               strand2_ID = crossover['second_strand_ID']
               strand2_bindex = crossover['second_strand_base_index']
               print("    %d: vhelix_base_index: %d  strand1_ID: %d  strand1_bindex: %d" %  
                   (n, vhelix_base_index, strand1_ID, strand1_bindex))
               if strand2_ID:
                   print("       strand2_ID: %d  strand2_bindex: %d" % (strand2_ID, strand2_bindex))
               n += 1
           #_for crossover in crossovers
        #__for connection in helix_connectivity
    #__for vhelix in vhelix_list

    print("")
    print("========================= strands =========================") 
    strand_list = json_data[VisJsonFields.STRANDS]
    print(">>> number of strands %d" % len(strand_list)) 
    domains_list = json_data[VisJsonFields.DOMAINS]
    for strand in strand_list:
        id = strand["id"]
        is_scaffold = strand[VisJsonFields.IS_SCAFFLOD]
        bases = strand[VisJsonFields.BASES]
        vhelix_list = strand[VisJsonFields.VIRTUAL_HELICES]
        strand_domain_list = strand[VisJsonFields.DOMAINS]
        num_strand_domain_bases = 0
        for domain_id in strand_domain_list:
            domain = domains_list[domain_id]
            num_bases = domain["number_of_bases"]
            num_strand_domain_bases += num_bases 
        color = strand[VisJsonFields.COLOR]
        print("")
        print(">>> strand %s: scaffold: %s nbases: %d nvhelix: %d " % (id, is_scaffold, len(bases), len(vhelix_list)))
        print("              number of domains: %d  domains: %s" % (len(strand_domain_list), str(strand_domain_list)))
        print("              color: %g %g %g" % (color[0],color[1],color[2])) 
        print("              bases: ") 
        print("              total number of domain bases: %d" % num_strand_domain_bases ) 
        n = 0
        for base in bases:
            id = base['id']
            coord = base['coordinates']
            seq = base['sequence']
            print("              %d: id %4d  coord (%6g, %6g, %6g) seq %s" % (n, id, coord[0],coord[1],coord[2], seq))
            n += 1

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

        print("")
        print(">>> domain %s: strand: %d  connected to strand: %d  connected to domain: %d " % (id, strand_id, connected_strand, 
            connected_domain))
        print("              nbases: %d  bases %s " % (len(base_list), str(base_list)))
        print("              start_pos: (%g %g %g) end_pos (%g %g %g) " % (start_pos[0], start_pos[1], start_pos[2], 
            end_pos[0], end_pos[1], end_pos[2]))
        print("              orientation: (%g %g %g) " % (orientation[0], orientation[1], orientation[2]))
        print("              start_base_index: %d  end_base_index: %d " % (start_base_index, end_base_index))
        print("              color: (%g, %g, %g) " % (color[0],color[1],color[2]))

if __name__=="__main__":
    main()
    

