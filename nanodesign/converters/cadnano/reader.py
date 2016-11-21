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

""" 
This module is used to read caDNAno DNA origami design JSON files. 
"""
import csv
import json
import logging
import re
import sys
from .common import CadnanoLatticeName,CadnanoLatticeType,CadnanoJsonFields
from .design import CadnanoDesign,CadnanoVirtualHelix,CadnanoBase
from ...data.sequence import DnaSequence

class CadnanoReader(object):
    """The CadnanoReader class."""
    instance_count = 1

    def __init__(self):
        self._logger = logging.getLogger(__name__)   
        CadnanoReader.instance_count += 1

    def read_json(self,file_name):
        """Read a caDNAno DNA origami design JSON file.

        Args:
            file_name (string): The name of a caDNAno DNA origami design JSON file to read.

        Returns:
            CadnanoDesign: A CadnanoDesign object containing the information parsed from the input caDNAno 
                DNA origami design JSON file.
        """
        # read the json data from the file.
        # make sure to expand the path first so we can use relative paths or ~ expansion
        import os.path
        file_name = os.path.expanduser( file_name )
        self._logger.info("Reading caDNAno design file {}".format(file_name))
        with open(file_name) as json_file:
            json_data = json.load(json_file)

        # parse the json data into a CadnanoDesign object.
        design = self.parse_json_data(json_data)

        # Calculate the possible crossovers.
        design.calculate_possible_crossovers()

        return design 

    def read_csv(self,file_name):
        """Read a caDNAno DNA origami design sequence CSV file.

        Args:
            file_name (string): The name of a caDNAno DNA origami design CSV file to read.

        Returns:
            list (DnaSequence): A list of DnaSequence objects representing the sequence for staple or scaffold strands.
        """

        self._logger.info("Reading caDNAno design CSV file %s " % file_name)
        sequences = []
        num_lines = 1
        seq_start = [0,0]
        seq_end = [0,0]

        with open(file_name, 'rU') as csv_file:
            reader = csv.reader(csv_file, dialect='excel')
            for row in reader:
                if (row[0] == 'Start'):
                    continue
                #print( ">>> row=%s" % str(row))
                start = row[0]
                result = re.search(r"([0-9]+)\[([0-9]+)\]", start)
                seq_start[0] = int(result.group(1))
                seq_start[1] = int(result.group(2))
                end = row[1]
                result = re.search(r"([0-9]+)\[([0-9]+)\]", end)
                seq_end[0] = int(result.group(1))
                seq_end[1] = int(result.group(2))
                seq_letters = row[2]
                seq_len = row[3]
                #print( ">>> seq_start=%d %d" % (seq_start[0], seq_start[1]))
                seq = DnaSequence(seq_start, seq_end, seq_letters, seq_len)
                sequences.append(seq)
                num_lines += 1
        #__with open(file_name, 'rU') __
        self._logger.info("Number of sequences read %d" % len(sequences))
        return sequences

    def parse_json_data(self, json_data):
        """Parse caDNAno DNA origami design JSON data.

        Args:
            json_data: The data read from a caDNAno DNA origami design JSON file.
        """
        design = CadnanoDesign()
        num_bases = len(json_data[CadnanoJsonFields.VSTRANDS][0][CadnanoJsonFields.SCAF])
        self._logger.info("Number of bases in a virtual helix %d " % num_bases)
        design.max_base_id = num_bases-1 

        # determine lattice type
        if ( (num_bases % 21 == 0) and (num_bases % 32) == 0):
            lattice_type = CadnanoLatticeType.honeycomb
        elif (num_bases % 32 == 0):
            lattice_type = CadnanoLatticeType.square
        elif (num_bases % 21) == 0:
            lattice_type = CadnanoLatticeType.honeycomb
        else:
            lattice_type = CadnanoLatticeType.honeycomb
        self._logger.info("Lattice type %s " % CadnanoLatticeType.names[lattice_type])
        design.lattice_type = lattice_type

        # parse helix information 
        num_scaffold_bases = 0;
        max_row_json = max_col_json = 0
        num_helix = 0
        for json_helix in json_data[CadnanoJsonFields.VSTRANDS]:
            num = int(json_helix[CadnanoJsonFields.NUM])
            row = int(json_helix[CadnanoJsonFields.ROW])
            col = int(json_helix[CadnanoJsonFields.COL])
            deletions = json_helix[CadnanoJsonFields.SKIP]
            insertions = json_helix[CadnanoJsonFields.LOOP]
            helix = CadnanoVirtualHelix(num_helix, num, row, col, insertions, deletions)
            self._logger.debug("==================== virtual helix %d ==================== " % num_helix) 
            self._logger.debug("num %d " % num) 
            self._logger.debug("row %d " % row) 
            self._logger.debug("col %d " % col) 
            self._logger.debug("Number of insertions %d " % insertions.count(-1))
            self._logger.debug("Number of deletions %d " % deletions.count(-1))

            # Check for a vhelix with no bases.
            scaffold = json_helix[CadnanoJsonFields.SCAF]
            staples = json_helix[CadnanoJsonFields.STAP]
            if scaffold.count([-1,-1,-1,-1]) + staples.count([-1,-1,-1,-1]) == len(scaffold) + len(staples): 
                continue

            # Parse scaffold information
            num_scaffold = len(scaffold)
            for json_base in scaffold:
                initial_strand = json_base[0]
                initial_base   = json_base[1]
                final_strand   = json_base[2]
                final_base     = json_base[3]
                self._logger.debug("scaffold base=" + str(json_base))
                base = CadnanoBase(initial_strand, initial_base, final_strand, final_base)
                helix.scaffold_strands.append(base)
                bsum = sum(json_base)
                if ( bsum != -4 ):
                    num_scaffold_bases += 1;
            #__for json_base in scaffold

            # parse staple information 
            num_staples = len(staples)
            for json_base in staples:
                initial_strand = json_base[0]
                initial_base   = json_base[1]
                final_strand   = json_base[2]
                final_base     = json_base[3]
                self._logger.debug("staple base=" + str(json_base))
                base = CadnanoBase(initial_strand, initial_base, final_strand, final_base)
                helix.staple_strands.append(base)
                bsum = sum(json_base)
            #__for json_base in staples

            design.helices.append(helix)

            self._logger.debug("Number of scaffold bases %d " % num_scaffold_bases) 
            self._logger.debug("Number of deletions %d " % deletions.count(-1))

            # parse staple colors
            helix.staple_colors = json_helix[CadnanoJsonFields.STAP_COLORS]

            max_row_json = max(max_row_json, int(json_helix[CadnanoJsonFields.ROW])+1)
            max_col_json = max(max_col_json, int(json_helix[CadnanoJsonFields.COL])+1)
            num_helix += 1
        #__for json_helix 

        design.max_row = max_row_json
        design.max_col = max_col_json
        self._logger.info("Number of virtual helices read %d " % num_helix) 
        self._logger.info("Number of scaffold bases %d " % num_scaffold_bases) 
        self._logger.info("Maximum number of lattice rows %d " % max_row_json) 
        self._logger.info("Maximum number of lattice columns %d " % max_col_json) 
        return design

