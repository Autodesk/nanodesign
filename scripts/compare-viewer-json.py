#!/usr/bin/python
# -*- coding: utf-8 -*-

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

"""This script is used to compare two Autodesk Nanodesign viewer JSON files.

   An Autodesk Nanodesign viewer JSON file contains data defining a DNA
   structure: bases, helices, strands, and domains. It is used by the Autodesk
   Nanodesign viewer to visualize a DNA structure translated from a DNA
   structure design file (e.g. caDNAno design file).

   Classes are defined to store data for the components (bases, domains, etc.)
   of the DNA structure. Each class defines methods to parse and compare
   component JSON data. The JSON field names specific to each component are
   stored in its class and used to parse its JSON data upon object
   initialization.

"""

import json
import os
import sys
import logging
import math

class Domain(object):
    """ The Domain class defines the data and methods for a DNA structure domain.
    """
    # JSON field names for domain data.
    BASES = 'bases'
    COLOR = 'color'
    CONN_DOM = 'connected_domain'
    CONN_STRAND = 'connected_strand'
    END_BASE_INDEX = 'end_base_index'
    END_POS = 'end_position'
    ID = 'id'
    MELTING_TEMP = 'melting_temperature'
    NUM_BASES = 'number_of_bases'
    ORIENTATION = 'orientation'
    START_BASE_INDEX = 'start_base_index'
    START_POS = 'start_position'
    STRAND_ID = 'strand_id'

    def __init__(self, json_data):
        """ Create a Domain object from JSON data. """
        self.bad_data = False
        self.num_errors = 0
        # Extract JSON data.
        try:
            self.id = json_data[self.ID]
            self.number_of_bases = json_data[self.NUM_BASES]
            self.base_list = json_data[self.BASES]
            self.strand_id = json_data[self.STRAND_ID]
            self.connected_strand = json_data[self.CONN_STRAND]
            self.connected_domain = json_data[self.CONN_DOM]
            self.start_base_index = json_data[self.START_BASE_INDEX]
            self.end_base_index = json_data[self.END_BASE_INDEX]
            self.orientation = json_data[self.ORIENTATION]
            self.start_pos = json_data[self.START_POS]
            self.end_pos = json_data[self.END_POS]
            self.color = json_data[self.COLOR]
            self.melting_temperature = json_data[self.MELTING_TEMP]
        except KeyError as e:
            logger.error("Domain data is missing the JSON %s field: %s" % (str(e), str(json_data)))
            self.bad_data = True

    def compare(self, domain):
        """ Compare this domain object with the given domain object. """
        if self.bad_data or domain.bad_data:
            return
        logger.info("Compare domain id %d " % self.id)
        self.num_errors += compare_dicts(self.__dict__, domain.__dict__)

        # Compare base ID lists.
        if not compare_lists(self.BASES, self.base_list, domain.base_list, check_order=True):
            self.num_errors += 1

        # Compare float lists.
        if not compare_float_lists(self.orientation,domain.orientation):
           logger.error("'%s' values are not equal: %s != %s" % (self.ORIENTATION, self.orientation, domain.orientation))
           self.num_errors += 1
        if not compare_float_lists(self.start_pos,domain.start_pos):
           logger.error("'%s' values are not equal: %s != %s" % (self.START_POS, self.start_pos, domain.start_pos))
           self.num_errors += 1
        if not compare_float_lists(self.end_pos,domain.end_pos):
           logger.error("'%s' values are not equal: %s != %s" % (self.END_POS, self.end_pos, domain.end_pos))
           self.num_errors += 1

class Base(object):
    """ The Base class defines the data and methods for a DNA base.
    """
    # JSON field names for base data.
    COORDS = 'coordinates'
    ID = 'id'
    SEQ = 'sequence'
    HELIX = 'h'
    POSITION = 'p'

    def __init__(self, json_data, id=0):
        """ Create a Base object from JSON data. """
        self.id = id
        self.bad_data = False
        self.num_errors = 0 
        try:
            # self.id = json_data[self.ID]    don't compare IDs, they can be different.
            self.coord = json_data[self.COORDS]
            self.sequence = json_data[self.SEQ]
            self.h = json_data[self.HELIX]
            self.p = json_data[self.POSITION]
        except KeyError as e:
            logger.error("Base data is missing the JSON %s field: %s" % (str(e), str(json_data)))
            self.bad_data = True

    def compare(self, base):
        """ Compare this base object with the given base object. """
        if self.bad_data or base.bad_data:
            return
        self.num_errors += compare_dicts(self.__dict__, base.__dict__)
        if not compare_float_lists(self.coord,base.coord):
           logger.error("Base number %d '%s' values are not equal: %s != %s" % (self.id, self.COORDS, self.coord, base.coord))
           self.num_errors += 1

class Strand(object):
    """ The Strand class defines the data and methods for a DNA strand. """
    # JSON field names for strand data.
    BASES = 'bases'
    DOMAINS = 'domains'
    ID = 'id'
    IS_SCAFFLOD = 'is_scaffold'
    IS_CIRCULAR = 'is_circular'
    NUM_BASES = 'number_of_bases'
    VIRTUAL_HELICES = 'virtual_helices'

    def __init__(self, json_data):
        """ Create a Strand object from JSON data. """
        self.bad_data = False
        self.num_errors = 0
        try:
            self.id = json_data[self.ID]
            self.is_scaffold = json_data[self.IS_SCAFFLOD]
            self.is_circular = json_data[self.IS_CIRCULAR]
            self.number_of_bases = json_data[self.NUM_BASES]
            self.bases = []
            for i,base in enumerate(json_data[self.BASES]):
                self.bases.append( Base(base,i) )
            self.helix_list = [] 
            for id in json_data[self.VIRTUAL_HELICES]:
                self.helix_list.append(id)
            self.domain_list = [] 
            for id in json_data[self.DOMAINS]:
                self.domain_list.append(id)
        except KeyError as e:
            logger.error("Strand data is missing the JSON %s field: %s" % (str(e), str(json_data)))
            self.bad_data = True

    def compare(self, strand):
        """ Compare this strand object with the given strand object. 

            Arguments:
                strand (DnaStrand): The strand object to compare.

            Strands are compared using
                1) The number of bases in the strands.
                2) The list of domains IDs. This may change because the integer domain 
                   IDs does not really matter.
                3) The bases in the strand. 
                4) The list of helix IDs. 
        """
        if self.bad_data or strand.bad_data:
            return
        logger.info("Compare strand id %d " % self.id)

        # Compare attributes that are not a list or a dict.
        self.num_errors += compare_dicts(self.__dict__, strand.__dict__)

        # Compare the list of domain IDs. 
        if not compare_lists(self.DOMAINS, self.domain_list, strand.domain_list, check_order=True):
            self.num_errors += 1

        # Compare the list of helix IDs. 
        if not compare_lists(self.VIRTUAL_HELICES, self.helix_list, strand.helix_list, check_order=True):
            self.num_errors += 1

        # Check length of bases.
        len1 = len(self.bases)
        len2 = len(strand.bases)
        if len1 != len2: 
            logger.error("'%s' the number of bases are not equal: %d != %d" % (self.BASES, len1, len2))
            self.num_errors += 1
            return

        # Check strand base data. 
        base_map1 = {}
        for base in self.bases:
            id = (base.h,base.p)
            base_map1[id] = base 
        base_map2 = {}
        for base in strand.bases:
            id = (base.h,base.p)
            base_map2[id] = base 
        id_list = []
        for id in base_map1:
            if id not in base_map2:
                id_list.append(id)
        if id_list:
            logger.error("'%s' IDs are different: %s not found in second file." % (self.BASES, id_list))
            self.num_errors += 1
            return 

        # Check base data.
        for id, base in base_map1.items():
            base.compare(base_map2[id])
            self.num_errors += base.num_errors

class HelixCrossover(object):
    """ The HelixCrossover class defines the data and methods for a helix crossover.
    """
    # JSON field names for helix crossover data.
    BASE_INDEX = 'vhelix_base_index'
    FIRST_STRAND_ID = 'first_strand_ID'
    FIRST_STRAND_BASE_INDEX = 'first_strand_base_index'
    SECOND_STRAND_ID = 'second_strand_ID'
    SECOND_STRAND_BASE_INDEX = 'second_strand_base_index'

    def __init__(self, json_data):
        """ Create a HelixCrossover object from JSON data. """
        self.bad_data = False
        self.num_errors = 0
        try:
             self.vhelix_base_index = json_data[self.BASE_INDEX]
             self.first_strand_ID = json_data[self.FIRST_STRAND_ID]
             self.first_strand_base_index = json_data[self.FIRST_STRAND_BASE_INDEX]
             self.second_strand_ID = json_data[self.SECOND_STRAND_ID]
             self.second_strand_base_index = json_data[self.SECOND_STRAND_BASE_INDEX]
        except KeyError as e:
            logger.error("Helix crossover data is missing the JSON %s field: %s" % (str(e), str(json_data)))
            self.bad_data = True

    def compare(self, crossover):
        """ Compare this crossover object with the given crossover object. """
        if self.bad_data or crossover.bad_data: 
            return
        # Compare attributes that are not a list or a dict.
        self.num_errors += compare_dicts(self.__dict__, crossover.__dict__)

class HelixConnectivity(object):
    """ The HelixConnectivity class defines the data and methods for helix connectivity.
    """
    # JSON field names for helix connectivity data.
    ANGLE = 'angle'
    CROSSOVERS = 'crossovers'
    DIRECTION = 'direction'
    HELIX_ID = 'helix_id'
    HELIX_NUM = 'helix_num'

    def __init__(self, json_data):
        """ Create a HelixConnectivity object from JSON data. """
        self.num_errors = 0 
        self.bad_data = False
        try:
            self.helix_id = json_data[self.HELIX_ID]
            self.helix_num = json_data[self.HELIX_NUM]
            self.angle = json_data[self.ANGLE]
            self.dir = json_data[self.DIRECTION]
            self.crossovers = []
            crossovers = json_data[self.CROSSOVERS]
            if crossovers:
                for crossover in crossovers:
                    self.crossovers.append( HelixCrossover(crossover) )
                #_for crossover in crossovers
        except KeyError as e:
            logger.error("Helix connectivity data is missing the JSON %s field: %s" % (str(e), str(json_data)))
            self.bad_data = True

    def compare(self, helix_conn):
        """ Compare this helix connectivity object with the given helix connectivty object. """
        if self.bad_data or helix_conn.bad_data: 
            return
        # Compare attributes that are not a list or a dict.
        self.num_errors += compare_dicts(self.__dict__, helix_conn.__dict__)

        # Compare directions. 
        if not compare_float_lists(self.dir,helix_conn.dir):
           logger.error("'%s' values are not equal: %s != %s" % (self.DIRECTION, self.dir, helix_conn.dir))
           self.num_errors += 1

        # Compare the number of crossovers. 
        len1 = len(self.crossovers)
        len2 = len(helix_conn.crossovers)
        if len1 != len2: 
            logger.error("'%s' the number of crossovers are not equal: %d != %d" % (self.CROSSOVERS, len1, len2))
            self.num_errors += 1
            return

        # Check that crossovers have the same IDs.
        map1 = {}
        for crossover in self.crossovers:
            map1[crossover.vhelix_base_index] = crossover
        map2 = {}
        for crossover in helix_conn.crossovers:
            map2[crossover.vhelix_base_index] = crossover
        id_list = []
        for crossover in self.crossovers:
            if crossover.vhelix_base_index not in map2: 
                id_list.append(crossover.vhelix_base_index)
        if id_list:
            logger.error("'%s' IDs are different: %s not found in second file." % (self.CROSSOVERS, id_list))
            return 

        # Compare crossover data.
        for id,crossover in map1.items():
            logger.info("'%s' %s %d " % (self.CROSSOVERS, crossover.BASE_INDEX, id)) 
            crossover.compare(map2[id])
            self.num_errors += crossover.num_errors 

class Helix(object):
    """ The Helix class defines the data and methods for a helix.
    """
    # JSON field names for helix data.
    ID = 'id'
    BP_RISE = 'base_pair_rise'
    CADNANO_INFO = 'cadnano_info'
    CADNANO_HELIX_NUM = 'num'
    CADNANO_COL = 'col'
    CADNANO_ROW = 'row'
    DOMAINS = 'domains'
    END_POS = 'end_position'
    HELIX_CONNECTIVITY = 'helix_connectivity'
    HELIX_DISTANCE = 'helix_distance'
    NUM_POSSIBLE_STAPLE_XOVERS = 'num_possible_staple_crossovers'
    NUM_POSSIBLE_SCAFFOLD_XOVERS = 'num_possible_scaffold_crossovers'
    ORIENTATION = 'orientation'
    START_POS = 'start_position'

    def __init__(self, json_data):
        """ Create a Helix object from JSON data."""
        self.num_errors = 0
        self.bad_data = False
        try:
            self.id = json_data[self.ID]
            cadnano_info = json_data[self.CADNANO_INFO]
            self.num = cadnano_info[self.CADNANO_HELIX_NUM]
            self.row = cadnano_info[self.CADNANO_ROW]
            self.col = cadnano_info[self.CADNANO_COL]
            self.helix_distance = json_data[self.HELIX_DISTANCE]
            self.base_pair_rise = json_data[self.BP_RISE]
            self.orientation = json_data[self.ORIENTATION]
            self.end_pos = json_data[self.END_POS]
            self.start_pos = json_data[self.START_POS]
            # Set domain list.
            self.domain_ids = []
            domains = json_data[self.DOMAINS]
            for domain_id in domains:
                self.domain_ids.append(domain_id)
            #__for domain_id in domains
            self.num_possible_staple_crossovers = json_data[self.NUM_POSSIBLE_STAPLE_XOVERS]
            self.num_possible_scaffold_crossovers = json_data[self.NUM_POSSIBLE_SCAFFOLD_XOVERS]
            # Set helix connectivity.
            self.helix_connectivity = []
            for connection in json_data[self.HELIX_CONNECTIVITY]:
                if not connection:
                    continue
                self.helix_connectivity.append( HelixConnectivity(connection) )
        except KeyError as e:
            logger.error("Helix data is missing the JSON %s field: %s" % (str(e), str(json_data)))
            self.bad_data = True

    def compare(self, helix):
        """ Compare this helix object with the given helix object. """
        if self.bad_data or  helix.bad_data:  
            return
        logger.info("Compare helix id %d  num %d" % (self.id, self.num))

        # Compare attributes that are not a list or a dict.
        self.num_errors += compare_dicts(self.__dict__, helix.__dict__)

        # Compare orientations, start and end positions. 
        if not compare_float_lists(self.orientation,helix.orientation):
            logger.error("'%s' values are not equal: %s != %s" % (self.ORIENTATION, self.orientation, helix.orientation))
            self.num_errors += 1
        if not compare_float_lists(self.start_pos,helix.start_pos):
            logger.error("'%s' values are not equal: %s != %s" % (self.START_POS, self.start_pos, helix.start_pos))
            self.num_errors += 1
        if not compare_float_lists(self.end_pos,helix.end_pos):
            logger.error("'%s' values are not equal: %s != %s" % (self.END_POS, self.end_pos, helix.end_pos))
            self.num_errors += 1

        # Compare domain ID lists.
        if not compare_lists("domains", self.domain_ids, helix.domain_ids, check_order=True):
            self.num_errors += 1

        # Compare helix connectivity. 
        self.compare_hconn(helix)

    def compare_hconn(self, helix):
        """ Compare helix connectivity. """
        len1 = len(self.helix_connectivity) 
        len2 = len(helix.helix_connectivity)
        if len1 != len2: 
            logger.error("'%s' the number of connected helices are not equal: %d != %d" % (self.HELIX_CONNECTIVITY, len1, len2))
            self.num_errors += 1
            return

        # Create a helix map indexed using (ID,num) tuple.
        conn_map1 = {}
        for connection in self.helix_connectivity:
            cindex = (connection.helix_id,connection.helix_num)
            conn_map1[cindex] = connection 
        conn_map2 = {}
        for connection in helix.helix_connectivity:
            cindex = (connection.helix_id,connection.helix_num)
            conn_map2[cindex] = connection 

        # Check IDs and nums.
        id_list = []
        for cindex in conn_map1:
            if cindex not in conn_map2:
                id_list.append(cindex)
        if id_list:
            logger.error("'%s' ID/num pairs are different: %s not found in second file." % (self.HELIX_CONNECTIVITY, id_list))
            self.num_errors += 1
            return 

        # Check connectivity data.
        for cindex,conn in conn_map1.items():
            logger.info("'%s' id %d  num %d " % (self.HELIX_CONNECTIVITY, cindex[0], cindex[1])) 
            conn.compare(conn_map2[cindex])
            self.num_errors += conn.num_errors

class DnaStructure(object):
    """ The DnaStructure class defines the data and methods for a DNA structure."""
    # JSON field names for DNA structure data.
    DOMAINS = "domains"
    LATTICE_TYPE = "lattice_type"
    STRANDS = "strands"
    VIRTUAL_HELICES = "virtual_helices"

    def __init__(self, name):
        self.name = name
        self.lattice_type = None 
        self.helix_list = []
        self.strand_list = []
        self.domain_list = []
        self.num_errors = 0

    def read_json(self, file_name):
        """ Read a JSON file. """
        logger.info("Read file %s" % file_name)
        with open(file_name) as file:
            json_data = json.load(file)
        self.parse_json(json_data)

    def parse_json(self, json_data):
        """ Parse the JSON data for a DNA structure. """
        if self.LATTICE_TYPE not in json_data:
            logger.error("The '%s' field is missing from the JSON file." % self.LATTICE_TYPE)
        else:
            self.lattice_type = json_data[self.LATTICE_TYPE]
            logger.info("Lattice type %s" % self.lattice_type)
        self.parse_helix_data(json_data)
        self.parse_strand_data(json_data)
        self.parse_domain_data(json_data)

    def parse_helix_data(self, json_data):
        """ Parse the JSON data for helices. """
        if self.VIRTUAL_HELICES not in json_data:
            logger.error("The '%s' field is missing from the JSON file." % self.VIRTUAL_HELICES)
            return 
        helix_list = json_data[self.VIRTUAL_HELICES]
        logger.info("Number of virtual helices %d" % len(helix_list)) 
        for helix_data in helix_list:
            self.helix_list.append( Helix(helix_data) )

    def parse_strand_data(self, json_data):
        """ Parse the JSON data for strands. """
        if self.STRANDS not in json_data:
            logger.error("The '%s' field is missing from the JSON file." % self.STRANDS)
            return 
        strand_list = json_data[self.STRANDS]
        logger.info("Number of strands %d" % len(strand_list))
        for strand in strand_list:
            self.strand_list.append( Strand(strand) )

    def parse_domain_data(self, json_data):
        """ Parse the JSON data for domains. """
        if self.DOMAINS not in json_data:
            logger.error("The '%s' field is missing from the JSON file." % self.DOMAINS)
            return 
        domains_list = json_data[self.DOMAINS]
        logger.info("Number of domains %d" % len(domains_list)) 
        for domain in domains_list:
            self.domain_list.append( Domain(domain) )

    def compare(self, dna_structure):
        """ Compare this object with the given DnaStructure object."""
        logger.info("")
        logger.info("Compare %s to %s:" % (self.name, dna_structure.name)) 
        if self.lattice_type != dna_structure.lattice_type: 
            logger.error("Lattice types not equal: %s != %s" % (self.lattice_type, dna_structure.lattice_type))
            self.num_errors += 1
        self.compare_helix(dna_structure)
        self.compare_strand(dna_structure)
        self.compare_domain(dna_structure)
        if self.num_errors != 0:
            logger.error("Files are not identical: %d errors detected." % self.num_errors)
            return False
        else:
            logger.info("Files are identical.")
            return True

    def compare_domain(self, dna_structure):
        """ Compare domain data. """
        logger.info("==================== compare domain data ====================")
        if len(self.domain_list) != len(dna_structure.domain_list):
            logger.error("Number of domains are not equal: %d != %d" % (len(self.domain_list), len(dna_structure.domain_list)))
            self.num_errors += 1
        dmap1 = self.get_domain_map()
        dmap2 = dna_structure.get_domain_map()
        id_list = []
        for id in dmap1:
            if id not in dmap2:
                id_list.append(id)
        #__for id in dmap1
        if id_list:
            logger.error("Domain IDs are not equal: %s not found in second file domain IDs." % id_list)
            self.num_errors += 1
            return

        # Compare domain data.
        for id,domain in dmap1.items(): 
            domain.compare(dmap2[id])
            self.num_errors += domain.num_errors 

    def get_domain_map(self):
        """ Create a dictionary mapping domain IDs to domain objects. """
        dmap = {}
        for domain in self.domain_list:
            dmap[domain.id] = domain
        return dmap 

    def compare_strand(self, dna_structure):
        logger.info("==================== compare strand data ====================")
        if len(self.strand_list) !=  len(dna_structure.strand_list):
            logger.error("Number of strands are not equal: %d != %d" % (len(self.strand_list), len(dna_structure.strand_list)))
            self.num_errors += 1
        smap1 = self.get_strand_map()
        smap2 = dna_structure.get_strand_map()
        id_list = []
        for id in smap1: 
            if id not in smap2: 
                id_list.append(id)
        #__for id in smap1
        if id_list:
            logger.error("Strand IDs are not equal: %s not found in second file strand IDs." % id_list)
            self.num_errors += 1
            return

        # Compare strand data.
        for id,strand in smap1.items(): 
            strand.compare(smap2[id])
            self.num_errors += strand.num_errors 

    def get_strand_map(self):
        """ Create a dictionary mapping strand IDs to strand objects. """
        smap = {}
        for strand in self.strand_list:
            smap[strand.id] = strand
        return smap 

    def compare_helix(self, dna_structure):
        logger.info("==================== compare helix data ====================")
        if len(self.helix_list) !=  len(dna_structure.helix_list):
            logger.error("Number of helices are not equal: %d != %d" % (len(self.helix_list), len(dna_structure.helix_list)))
            self.num_errors += 1

        # Check that the helix 'num' fields match.
        map1 = self.get_helix_map()
        map2 = dna_structure.get_helix_map()
        num_list = []
        for num in map1.keys():
            if num not in map2:
                num_list.append(id)
        #__for num in map1.keys()
        if num_list:
            logger.error("Helix nums are not equal: %s not found in second file helix nums." % num_list)
            self.num_errors += 1
            return
 
        # Compare helix data.
        for num, helix in map1.items():
            helix.compare(map2[num])
            self.num_errors += helix.num_errors 

    def get_helix_map(self):
        """ Create a dictionary mapping helix num to helix objects. """
        hmap = {}
        for helix in self.helix_list:
            hmap[helix.num] = helix
        return hmap 

def compare_dicts(dict1, dict2):
    """ Compare the non-list, non-dict entries of two dictionaries. """
    num_errors = 0
    for attr, value1 in dict1.iteritems():
        value2 = dict2[attr]
        if (type(value1) in (bool,int,float,str,unicode)) and (value1 != value2):
            logger.error("'%s' values are not equal: %s != %s" % (attr, str(value1), str(value2)))
            num_errors += 1
    return num_errors

def compare_float_lists(list1, list2, tol=1e-9):
    """ Compare two float lists. 
        This function is used to compare two float lists, usually representing a vector. 
    """
    if len(list1) != len(list2):
        return False

    for val1,val2 in zip(list1, list2):
        if abs(val1-val2) > tol:
            return False
    #__for val1,val2 in zip(list1, list2)

    return True

def compare_lists(name, list1, list2, check_order):
    """ Compare two integer lists. 

        This function is used to compare two integer lists, usually representing a list of
        component IDs (e.g. domain IDs). The order of the IDs may or may not be important.
    """
    if len(list1) != len(list2):
        logger.error("'%s' lists are not the same size: %d != %d" % (name, len(list1), len(list2)))
        return False

    # Check that the lists contain the same elements.
    if (False):
        ldiff = set(list1).symmetric_difference(set(list2))
        if ldiff:
            logger.error("'%s' lists contain different values: %s" % (name, list(ldiff)))
            return False

        # Check if the ordering of the lists are the same.
        if check_order and cmp(list1,list2):
            logger.error("'%s' list values are not in the same order: %s != %s" % (name, list1, list2)) 
            return False

    return True

def main():
    global logger
    logger = logging.getLogger('compare-viewer-json')

    if (len(sys.argv) != 3):
        sys.stderr.write("**** ERROR: Wrong number of arguments.\n") 
        sys.stderr.write("Usage: compare-viewer-json.py <file1> <file2> \n")
        sys.exit(1)

    # Read the 1st file and create a DNA structure object from it.
    file_name1 = sys.argv[1]
    dna_structure1 = DnaStructure(file_name1)
    dna_structure1.read_json(file_name1)
    logger.info("")

    # Read the 2nd file and create a DNA structure object from it.
    file_name2 = sys.argv[2]
    dna_structure2 = DnaStructure(file_name2)
    dna_structure2.read_json(file_name2)

    # Compare the data.
    if not dna_structure1.compare(dna_structure2):
        sys.exit(1)

if __name__ == "__main__":
    main()
