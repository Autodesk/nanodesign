#!/usr/bin/env python
""" 
This module is used to create a DNA structure from a caDNAno DNA origami design file. 

The topology table provides connectivity information between all the bases in a model. The coordinates each base on the DNA
helix is also computed and stored in the table. The positions of base pairs along the DNA helix axes with their orientations are 
also computed and stored in separate data structures. 

This code is based on a direct translation of the set MATLAB scripts to convert a caDNAno design to a CanDo 
.cndo file from Mark Bathe's Laboratory for Computational Biology & Biophysics at MIT.
"""
import os
import sys
import csv
import logging
import re
import string
import sys
import time
import numpy as np
from math import sqrt,cos,acos,sin,asin,pi
from dna_sequence_data import dna_sequence_data

from design import CadnanoDesign,CadnanoVirtualHelix,CadnanoBase
from reader import CadnanoReader 
from common import CadnanoLatticeType
from nanodesign.base import DnaBase 
from nanodesign.strand import DnaStrand
from nanodesign.dna_structure import DnaStructure,DnaStructureHelix

# The number of columns in the intermediate structure topology table.
TOPOLOGICAL_TABLE_NUM_COLUMNS = 18

class StrandType:
    SCAFFOLD = 0
    STAPLE   = 1

class Tindex:
    """ The indexes into the structure topology table. """
    ID = 0
    CURRENT_BASE   = 1
    FIVE_NEIGHBOR  = 4
    THREE_NEIGHBOR = 7
    WC_NEIGHBOR    = 10
    COORD          = 14
    DELETION       = 16
    INSERTION      = 17

class CadnanoConvertDesign(object):
    """ The CadnanoConvertDesign class creates a DNA structure from a caDNAno DNA origami design. 

        Attributes:
            dna_structure (DnaStructure): A DnaStructure object. 
    """

    def __init__(self):
        self._logging_level = logging.INFO
        self._setup_logging()
        self._timer = _Timer()
        self.dna_structure = None 

    def _set_logging_level(self,level):
        """Set logging level."""
        self._logger.setLevel(level)

    def _setup_logging(self):
        """ Set up logging."""
        self._logger = logging.getLogger('cadnano:convert_design')
        self._logger.setLevel(self._logging_level)

        # Create a console handler and set format.
        console_handler = logging.StreamHandler()
        #formatter = logging.Formatter('%(asctime)s [%(name)s] %(levelname)s - %(message)s')
        formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        self._logger.addHandler(console_handler)

    def create_structure(self,design):
        """ Create a DNA structure from a caDNAno DNA origami design.

            Arguments:
                design (CadnanoDesign): A caDNAno DNA origami design.

            Returns:
                DnaStructure: The DNA structure object containing topological and geometric information for the
                              caDNAno DNA origami design model.
        """

        self._timer.start()
        structure_topology, dnode, triad, id_nt_0, helices = self._create_structure_topology_and_geometry(design)
        num_bases = structure_topology.shape[0]
        self._logger.info("Number of bases in structure_toplogy: %d " % num_bases)
        self._logger.info("Time to create structure topology table: %s " % self._timer.finish())
        self._logger.info("Number of bases in helix axis: %d " % dnode.shape[0])

        # Create the nucleotide map table.
        self._timer.start()
        id_nt = self._create_nt_map_table(structure_topology, dnode, id_nt_0)
        self._logger.info("Time to create map table: %s " % self._timer.finish())

        # Renumber base IDs.
        for i in xrange(0,num_bases):
            structure_topology[i,0] = i+1

        # Create the DNA topology table.
        self._timer.start()
        dna_topology = self._create_dna_topology(structure_topology)
        self._logger.info("Time to create topology table: %s " % self._timer.finish())

        #if (True):
        if (False):
            print("---------- dna_topology before inserts/deletes ----------")
            for i in xrange(0,len(dna_topology)):
                base = dna_topology[i];
                print("%d:  id=%d  up=%d  down=%d  across=%d  h=%d  p=%d  isScaf=%d" % 
                    (i+1, base.id, base.up, base.down, base.across, base.h, base.p, base.is_scaf) )

        # store the topology and geometry information into a DnaStructure object.
        self.dna_structure = DnaStructure()
        self.dna_structure.base_connectivity = dna_topology
        self.dna_structure.helix_axis_nodes = dnode 
        self.dna_structure.helix_axis_frames = triad
        self.dna_structure.id_nt = id_nt
        self.dna_structure.lattice_type = design.lattice_type
        self.dna_structure.add_structure_helices(helices)

        # Remove deleted bases
        delete_and_insert = False
        if (delete_and_insert):
            dna_topology, dnode, triad, id_nt = self._delete_bases(dna_topology, dnode, triad, id_nt)

        # Add inserted bases
        if (delete_and_insert):
            dna_topology, dnode, triad, id_nt = self._insert_bases(dna_topology, dnode, triad, id_nt)

        # Generate strands 
        dna_topology,strands = self._build_strands(dna_topology)
        self.dna_structure.strands = strands
        self._logger.info("Number of strands: %d " % len(strands)) 

        # Calculate staple ends 
        staple_ends = self._calculate_staple_ends(dna_topology)
        self.dna_structure.staple_ends = staple_ends
 
        return self.dna_structure

    def _create_dna_topology(self,structure_topology):
        """ Create the DNA topology table. """
        num_bases = structure_topology.shape[0]
        dna_topology = []

        for i in xrange(0,num_bases):
            id = structure_topology[i][0]
            base = DnaBase(id)

            #----- get 5' neighbor -----#
            five_neigh = structure_topology[i,4:7]
            curr_bases = structure_topology[:,1:4]
            tmp = _find_row(five_neigh, curr_bases)

            if (len(tmp) == 1):
                base.up = int(tmp[0])+1
            elif (len(tmp) == 0):
                base.up = -1
            else:
                logger.error('Duplication.')

            #----- get 3' neighbor -----#
            three_neigh = structure_topology[i,7:10]
            curr_bases = structure_topology[:,1:4]
            tmp = _find_row(three_neigh, curr_bases)

            if (len(tmp) == 1):
                base.down = int(tmp[0])+1
            elif (len(tmp) == 0):
                base.down = -1
            else:
                logger.error('Duplication.')

            #----- Watson-Crick neighbor -----#
            wc_neigh = structure_topology[i,10:13]
            curr_bases = structure_topology[:,1:4]
            tmp = _find_row(wc_neigh, curr_bases)

            if (len(tmp) == 1):
                base.across = int(tmp[0])+1
            elif (len(tmp) == 0):
                base.across = -1
            else:
                self._logger.error('Duplication.')

            # helix, position, scaffold/staple
            base.h = structure_topology[i,1]
            base.p = structure_topology[i,2]

            if (structure_topology[i,3] == 0):
                base.is_scaf = True
            elif (structure_topology[i,3] == 1):
                base.is_scaf = False
            else:
                logger.error('Exception.')

            # coordinate
            # note: where is this set?
            #base.coord = structure_topology[i,13:16]

            # deletion
            base.skip = int(structure_topology[i,16])

            # insertion
            base.loop = int(structure_topology[i,17])

            dna_topology.append(base)
        #__for i in xrange(0,num_bases)__
        return dna_topology 

    def _delete_bases(self, dna_topology, dnode, triad, id_nt):
        """ Remove bases from dna_topology[]. deleted bases are given in the Base.skip field.
            Note that curr_skip is a zero-based index; the other base IDs in dna_topology are 1-based.
        """

        self._logger.debug("==================== delete bases ====================")
        num_bases = len(dna_topology)

        num_skip = 0
        for i in xrange(0,num_bases):
            if ( dna_topology[i].skip < 0):
               num_skip += 1

        self._logger.debug("Number of bases to delete: %d" % num_skip )
        if (num_skip == 0):
            return dna_topology, dnode, triad, id_nt 

        # create a list of indexes for skipped bases
        skip = np.zeros(num_bases, dtype=int)
        for i in xrange(0,num_bases):
            skip[i] = dna_topology[i].skip
        #self._logger.debug("skip: %s " % str(skip))

        # list of bases deleted
        deleted_bases = []
        deleted_bases_bp = []

        #----- iterated over skipped bases -----#
        nonzero_search = np.nonzero(skip)[0]
        while (nonzero_search.size != 0): 
            curr_skip = nonzero_search[0]
            #print(">>> curr_skip=%d" % curr_skip) 

            # find the four neighboring bases
            neighbor_up = dna_topology[curr_skip].up
            neighbor_down = dna_topology[curr_skip].down
            curr_skip_across = dna_topology[curr_skip].across

            if (curr_skip_across >= 0):
                neighbor_across_up = dna_topology[curr_skip_across-1].up
                neighbor_across_down = dna_topology[curr_skip_across-1].down
            else:
                neighbor_across_up = -1
                neighbor_across_down = -1

            # update base connectivity
            if (neighbor_up >= 0):
                dna_topology[neighbor_up-1].down = neighbor_down

            if (neighbor_down >= 0):
                dna_topology[neighbor_down-1].up = neighbor_up

            if (neighbor_across_up >= 0):
                dna_topology[neighbor_across_up-1].down = neighbor_across_down

            if (neighbor_across_down >= 0):
                dna_topology[neighbor_across_down-1].up = neighbor_across_up

            # update the list deleted bases
            if (curr_skip_across >= 0):
                deleted_bases.append(curr_skip)
                deleted_bases.append(curr_skip_across-1)
                ti = np.where(id_nt == curr_skip)[0]
                for i in xrange(0,ti.shape[0]):
                    deleted_bases_bp.append(ti[i])
                skip[curr_skip] = 0
                skip[curr_skip_across-1] = 0
            else:
                deleted_bases.append(curr_skip)
                skip[curr_skip] = 0
       
            # find the next base to be deleted
            nonzero_search = np.nonzero(skip)[0]
        #_while (nonzero_search.size != 0)_

        self._logger.debug("Deleted_bases: %s" % str(deleted_bases) )
        #self._logger.debug("Deleted_bases_bp: %s" % str(deleted_bases_bp) )

        # remove bases from topology
        for i in xrange(0,len(deleted_bases)):
            base = deleted_bases[i] - i
            #print(">>> delete %d %d %d %d" % ( dna_topology[base].id, dna_topology[base].up, dna_topology[base].down, 
            #      dna_topology[base].across))
            del dna_topology[base] 

        # remove bases from other arrays
        dnode = np.delete(dnode, deleted_bases_bp, 0)
        triad = np.delete(triad, deleted_bases_bp, 2)
        id_nt = np.delete(id_nt, deleted_bases_bp, 0)

        # update the base IDs
        dna_topology, id_nt = self._renumber_baseIDs(dna_topology, id_nt)
        self._logger.debug("Number of bases in topology table: %d " % num_bases)
        return dna_topology, dnode, triad, id_nt 

    def _insert_bases(self, dna_topology, dnode, triad, id_nt):
        """ Insert bases into dna_topology[]. Inserted bases are given in the Base.loop field. """
        dist_bp = 3.4       # Rise between two neighboring base-pairs (Angstrom)
        ang_bp = 360/10.5   # Twisting angle between two neighboring base-pairs (degree)

        self._logger.debug("==================== insert bases ====================")
        num_bases = len(dna_topology)
        dy = np.array([0, 0.1, 0], dtype=float)

        # these need to have same shape has slice of triad[]
        y_up_vec   = np.array([0,  1.0, 0], dtype=float)
        y_down_vec = np.array([0, -1.0, 0], dtype=float)

        # Check for inserts; if none then just return.
        num_insert = 0
        for i in xrange(0,num_bases):
            if (dna_topology[i].loop == 1):
                num_insert += 1
        self._logger.debug("Number of bases to insert: %d" % num_insert)
        if (num_insert == 0):
            return dna_topology, dnode, triad, id_nt 

        # create a loops array for searching
        loops = np.zeros((num_bases,), dtype=int);
        for i in xrange(0,num_bases):
            loops[i] = dna_topology[i].loop

        #----- iterated over inserted bases -----#
        num_added_inserts = 0 
        nonzero_search = np.nonzero(loops)[0]
        while (nonzero_search.size != 0):
            curr_loop = nonzero_search[0]
            #print("")
            #print("[insert_bases] ------------------------------------- " )
            #print("[insert_bases] curr_loop=%d" % curr_loop)

            # find the four neighboring bases
            curr_loop_across = dna_topology[curr_loop].across
            #print("[insert_bases] curr_loop_across=%d" % curr_loop_across)
            #print("[insert_bases] id_nt.shape=%s" % str(id_nt.shape))

            if (curr_loop_across >= 0):
                row_list,col_list = np.where(id_nt == curr_loop)
                ti = row_list[0]
                tj = col_list[0]
                #print("[insert_bases] ti=%s" % str(ti)) 
                #print("[insert_bases] tj=%s" % str(tj)) 
                neighbor_across_up = dna_topology[curr_loop_across-1].up
                # Make sure that the 3'-neighbor of the base 'currLoop' is to the right
                tmp_a = triad[:,2,ti]
                #print("[insert_bases] tmp_a=%s" % str(tmp_a)) 
                #print("[insert_bases] tmp_a.shape=%s" % str(tmp_a.shape)) 
                #print("[insert_bases] triad.shape=%s" % str(triad.shape)) 

                if (np.linalg.norm(tmp_a - y_up_vec) < 1e-10):  # reference axis e3 pointing to the right
                    #print("[insert_bases] ref axis points to right") 
                    if (tj == 1): # base 'currLoop' is a non-preferred base
                        tmp_nt = curr_loop
                        curr_loop = curr_loop_across-1
                        curr_loop_across = tmp_nt+1

                elif (np.linalg.norm(tmp_a - y_down_vec) < 1e-10): # reference axis e3 pointing to the left
                    #print("[insert_bases] ref axis points to left ") 
                    if (tj == 0): # base 'currLoop' is a preferred base
                        tmp_nt = curr_loop
                        curr_loop = curr_loop_across-1
                        curr_loop_across = tmp_nt+1
                else:
                    logger.error('Exception.')
            #__if (curr_loop_across >= 0)__

            # Find the four neighboring bases and insert them.
            neighbor_down = dna_topology[curr_loop].down
            if (curr_loop_across >= 0):
                neighbor_across_up = dna_topology[curr_loop_across-1].up
            else:
                neighbor_across_up = -1

            # insert the bases
            if (curr_loop_across >= 0):
                num_inserts = 2*loops[curr_loop]
            else:
                num_inserts = loops[curr_loop]

            #print("[insert_bases] num_inserts=%s" % str(num_inserts))
            #print("[insert_bases] neighbor_down=%d" % neighbor_down)
            #print("[insert_bases] curr_loop_across=%d" % curr_loop_across)
            #print("[insert_bases] neighbor_across_up=%d" % neighbor_across_up)

            #----- case I: insert dsDNA -----#

            if (curr_loop_across >= 0):
                # Create basepair information
                row_list,col_list = np.where(id_nt == curr_loop)
                ti = row_list[0]
                tj = col_list[0]
                #print("[insert_bases] ti=%s" % str(ti))
                #print("[insert_bases] tj=%s" % str(tj))

                # Positions and orientations
                dnode_1 = dnode[ti,:]
                #print("[insert_bases] dnode_1=%s" % str(dnode_1))
                #print("[insert_bases] dnode_1.shape=%s" % str(dnode_1.shape))

                triad_1 = triad[:,:,ti]
                triad_1 = triad_1.reshape(3,3)
                #print("[insert_bases] triad_1.shape=%s " % str(triad_1.shape))
                #print("[insert_bases] triad_1=%s" % str(triad_1))

                dnode_2 = dnode_1 + [0, dist_bp, 0]
                #print("[insert_bases] dnode_2=%s" % str(dnode_2))
                #print("[insert_bases] dnode_2.shape=%s " % str(dnode_2.shape))

                rot_mat = _vrrotvec2mat(y_up_vec, _deg2rad(ang_bp)) 
                #print("[insert_bases] rot_mat=%s " % str(rot_mat))
                #print("[insert_bases] rot_mat.shape=%s " % str(rot_mat.shape))
                triad_2 = np.dot(rot_mat,triad_1)

                #print("[insert_bases] triad_2=%s" % str(triad_2))
                #print("[insert_bases] triad_2.shape=%s " % str(triad_2.shape))

                [dnode_insert, triad_insert] = _bp_interp(dnode_1, triad_1, dnode_2, triad_2, num_inserts/2)

                #print("[insert_bases] dnode_insert=%s" % str(dnode_insert))
                #print("[insert_bases] triad_insert=%s" % str(triad_insert))

                # map table
                #id_nt_insert = reshape(nBase+(1:nInsert), 2, nInsert/2)';
                #id_nt_insert = np.arange(num_bases+1,num_bases+num_inserts+1)
                #id_nt_insert = np.zeros((num_inserts, 2), dtype=int)
                id_nt_insert = np.zeros((num_inserts/2, 2), dtype=int)
                for i in xrange(0,num_inserts/2):
                    #id_nt_insert[i][0] = num_bases + i+1 
                    #id_nt_insert[i][1] = num_bases + i+2 
                    id_nt_insert[i][0] = num_bases + i 
                    id_nt_insert[i][1] = num_bases + i+1 
                #print("[insert_bases] id_nt_insert=%s" % str(id_nt_insert))

                if (np.linalg.norm(triad_1[:,2] - y_up_vec) < 1e-10):
                    pass
                elif (np.linalg.norm(triad_1[:,2] - y_down_vec) < 1e-10):  # base 'currLoop' is non-preferred
                    #print("[insert_bases] base is non-preferred") 
                    #print("[insert_bases] id_nt_insert before flipls: %s " % str(id_nt_insert)) 
                    id_nt_insert = np.fliplr(id_nt_insert);
                    #print("[insert_bases] id_nt_insert after flipls: %s " % str(id_nt_insert)) 
                else:
                    logger.error('Exception.');

                # insert new basepairs
                dnode = np.concatenate((dnode, dnode_insert), axis=0);
                triad = np.concatenate((triad, triad_insert), axis=2);
                id_nt = np.concatenate((id_nt, id_nt_insert), axis=0);

                for i in xrange(0,num_inserts,2):
                    base1 = DnaBase(num_bases+i+1)
                    base2 = DnaBase(num_bases+i+2)
                    #print("[insert_bases] add two bases: id1=%d  id1=%d " % ( base1.id, base2.id) )
                    dna_topology.append(base1)
                    dna_topology.append(base2)
                    num_added_inserts += 2 

                    if (i == 0):
                        dna_topology[curr_loop].down = num_bases+1
                        dna_topology[curr_loop_across-1].up = num_bases+2
                        base1.up = curr_loop+1
                        base2.down = curr_loop_across
                    else:
                        base1.up = num_bases+i-2+1
                        base2.down = num_bases+i-1+1;

                    # [davep] is this right? -2 or -1?
                    if (i == num_inserts-2):
                        if (neighbor_down >= 0):
                            dna_topology[neighbor_down-1].up = num_bases + num_inserts-1

                        if (neighbor_across_up >= 0):
                            dna_topology[neighbor_across_up-1].down = num_bases + num_inserts

                        base1.down = neighbor_down
                        base2.up = neighbor_across_up
                    else:
                        base1.down = num_bases+i+2
                        base2.up = num_bases+i+3

                    base1.across = base2.id
                    base2.across = base1.id

                    base1.h = dna_topology[curr_loop].h
                    base2.h = dna_topology[curr_loop_across-1].h
                    base1.p = dna_topology[curr_loop].p
                    base2.p = dna_topology[curr_loop_across-1].p
                    base1.is_scaf = dna_topology[curr_loop].is_scaf
                    base2.is_scaf = dna_topology[curr_loop_across-1].is_scaf

                    # chang this Tue Feb  9 11:39:34 PST 2016
                    #base1.h = base2.h
                    #base2.h = base1.h
                    #base1.p = base2.p
                    #base2.p = base1.p
                    #base1.is_scaf = base2.is_scaf
                    #base2.is_scaf = base1.is_scaf

                    base1.coord = dna_topology[curr_loop].coord + dy*(i+1+1)/2
                    base2.coord = dna_topology[curr_loop_across-1].coord + dy*(i+1)/2
                    base1.skip = 0
                    base2.skip = 0
                    base1.loop = 0
                    base2.loop = 0
                #__for i in xrange(0,num_inserts)__

            # insert ssDNA

            else:
                for i in xrange(0,num_inserts):
                    base = DnaBase(num_bases+i+1)
                    dna_topology.append(base)
                    num_added_inserts += 1 
                    #print(">>> add one base: id=%d " % base.id )

                    if (i == 0):
                        dna_topology[curr_loop].down = num_bases+1
                        base.up = curr_loop+1
                    else:
                        base.up = num_bases+i-1+1

                    if (i == num_inserts-1):
                        dna_topology[neighbor_down-1].up = num_bases + num_inserts
                        base.down = neighbor_down
                    else:
                        base.down = num_bases+i+1+1

                    base.across = -1
                    base.skip = 0
                    base.loop = 0
                    base.h = dna_topology[curr_loop].h
                    base.p = dna_topology[curr_loop].p
                    base.is_scaf = dna_topology[curr_loop].is_scaf
                    base.coord = dna_topology[curr_loop].coord + dy*(i+1)
            #__if (curr_loop_across >= 0)__

            # update num_bases for creating new base IDs
            num_bases = len(dna_topology)

            # remove loop entries 
            if (curr_loop_across >= 0):
                loops[curr_loop] = 0
                loops[curr_loop_across-1] = 0
            else:
                loops[curr_loop] = 0

            # find the next base to be inserted 
            nonzero_search = np.nonzero(loops)[0]
        #_while (nonzero_search.size != 0)_

        #self._logger.debug("Number of bases inserted: %d " % num_added_inserts)
        print_on = True
        print_on = False
        if (print_on):
            num_bases = len(dna_topology)
            print("--------- new topology with inserts ---------")
            for i in xrange(0,num_bases):
                #print("%d: %d %d %d %d %f %f %f" % (i+1, dna_topology[i].id, dna_topology[i].up, 
                #      dna_topology[i].down, dna_topology[i].across, dna_topology[i].coord[0], 
                #      dna_topology[i].coord[1],  dna_topology[i].coord[2] )) 
                print("%d: %d %d %d %d" % (i+1, dna_topology[i].id, dna_topology[i].up, 
                      dna_topology[i].down, dna_topology[i].across ))

        self._logger.debug("Number of bases in topology table: %d " % num_bases)
        return dna_topology, dnode, triad, id_nt 

    def _renumber_baseIDs(self, dna_topology, id_nt):
        """ Renumber base IDs so that they are between 1 and len(dna_topology). This is needed when bases are added or deleted.
        """    
        num_bases = len(dna_topology)
        conn = np.zeros((num_bases,4), dtype=int)
        for i in xrange(0,num_bases):
            conn[i,:] = [dna_topology[i].id, dna_topology[i].up, dna_topology[i].down, dna_topology[i].across]

        # Sort the connectivty 
        sorted_indices = np.argsort(conn, axis=0)
        conn = conn[sorted_indices[:,0]]

        # Modify base IDs 
        # Note: be careful here. id_nt contains zero-indexed base IDs.
        id_nt_new = id_nt
        new_conn = -2*np.ones((num_bases,4), dtype=int)
        new_conn[conn == -1] = -1    # set entries in new_conn[] to -1 for -1 entries in conn[] 
        for i in xrange(0,num_bases):
            new_conn[conn == conn[i,0]] = i+1
            #id_nt_new[id_nt == conn[i,0]] = i+1
            id_nt_new[id_nt == conn[i,0]-1] = i

        # Create a new topology 
        new_dna_topology = []
        for i in xrange(0,num_bases):
            id = new_conn[i,0]
            base = DnaBase(id)
            base.up = new_conn[i,1]
            base.down = new_conn[i,2]
            base.across = new_conn[i,3]
            j = sorted_indices[i][0]
            base.h = dna_topology[j].h
            base.p = dna_topology[j].p
            base.is_scaf = dna_topology[j].is_scaf
            base.coord = dna_topology[j].coord
            base.skip = dna_topology[j].skip
            base.loop = dna_topology[j].loop
            new_dna_topology.append(base)
        #__for i in xrange(0,num_bases)
        return new_dna_topology, id_nt_new

    def _create_structure_topology_and_geometry(self,design):
        """ Create topological and geometrical information for a design. """
        structure_toplogy = np.empty(shape=(0,TOPOLOGICAL_TABLE_NUM_COLUMNS))
        dnode = np.empty(shape=(0,3),dtype=float)
        triad = np.empty(shape=(3,3,0),dtype=float)
        id_nt_0 = np.empty(shape=(0,6),dtype=float)
        num_vhelices = 0
        lattice_type = design.lattice_type
        row_list = []
        col_list = []
        structure_helices = [] 
        vhelices = design.helices
        self._logger.setLevel(logging.DEBUG)
        self._logger.debug("==================== create structure topology and geometry ====================")

        for vhelix in vhelices:
            self._logger.debug("---------- process virtual helix %d ----------" % num_vhelices);
            row = vhelix.row 
            col = vhelix.col 
            num = vhelix.num 
            row_list.append(row)
            col_list.append(col)
            self._logger.debug("num: %d " % num)
            self._logger.debug("row: %d " % row)
            self._logger.debug("col: %d " % col)
            structure_helix = DnaStructureHelix(num_vhelices)
            structure_helix.lattice_num = num
            structure_helix.lattice_row = row
            structure_helix.lattice_col = col

            if ( num % 2 == 0 ):
                structure_helix.scaffold_polarity = "5'"
                self._logger.debug("scaffold polarity 5' to 3'")
            else:
                structure_helix.scaffold_polarity = "3'"
                self._logger.debug("scaffold polarity 3' to 5'")

            # Set staple colors
            stap_colors = vhelix.staple_colors 
            for color in stap_colors:
                structure_helix.staple_colors.append(color)

            # Create data for a single helix
            helix_topology, dnode_0, triad_0, id, dnode_full = self._create_single_helix(vhelix, lattice_type)

            # Define the geometry of the helix by setting its end coordinates.
            structure_helix.end_coordinates[0] = dnode_0[0]
            structure_helix.end_coordinates[1] = dnode_0[-1]
            structure_helix.end_frames[:,:,0] = triad_0[:,:,0]
            structure_helix.end_frames[:,:,1] = triad_0[:,:,-1]
            structure_helix.helix_axis_nodes = dnode_full 

            # Append results to the global table
            structure_toplogy = np.concatenate((structure_toplogy, helix_topology), axis=0)
            dnode = np.concatenate((dnode, dnode_0), axis=0)
            triad = np.concatenate((triad, triad_0), axis=2)
            id_nt_0 = np.concatenate((id_nt_0, id), axis=0)

            structure_helices.append(structure_helix)
            num_vhelices += 1
        #__for vhelix in vhelices
        return structure_toplogy, dnode, triad, id_nt_0, structure_helices

    def _create_nt_map_table(self, structure_topology, dnode, id_nt_0):
        """ Create the map table (id_nt). """
        n_bp = dnode.shape[0]
        id_nt = np.zeros((n_bp, 2), dtype=int)
        for i in xrange(0,n_bp):
            # scaffold nucleotide
            tmp = _find_row(id_nt_0[i,0:3], structure_topology[:,1:4])
            id_nt[i,0] = tmp
            # staple nucleotide
            tmp = _find_row(id_nt_0[i,3:6], structure_topology[:,1:4])
            id_nt[i,1] = tmp;
        return id_nt 

    def _create_single_helix(self, vhelix, lattice_type):
        scaffolds = vhelix.scaffold_strands 
        staples = vhelix.staple_strands 
        deletions = vhelix.deletions 
        insertions = vhelix.insertions 
        row = vhelix.row 
        col = vhelix.col 
        num = vhelix.num 
        num_bases = len(scaffolds)
        #self._logger.debug("number of bases=%d" % num_bases )
        #self._logger.debug("vstrand.row=%d" % row)
        #self._logger.debug("vstrand.col=%d" % col)

        # generate the coordinates for the scaffold and staple bases
        dnode = np.empty((0, 3), dtype=float)
        triad = np.empty((3, 3,0), dtype=float)
        scaffold_coords, staple_coords, dnode_0, triad_0 = self._generate_coordinates(lattice_type, row, col, num, num_bases)
        # convert nm to Angstrom (scaffold_coords and staple_coords are in nm).
        #dnode_0 = dnode_0 * 10;

        # Add scaffold and strand information to the topology table. 
        uid = -1
        helix_topology = np.zeros((2*num_bases, TOPOLOGICAL_TABLE_NUM_COLUMNS), dtype=float)
        id_nt_0 = np.empty((0,6), dtype=float)

        for i in xrange(0,num_bases):
            current_scaffold = scaffolds[i] 
            current_staple = staples[i] 

            # if the base exists in the scaffold strand
            #if ((current_scaffold[0] >= 0) or (current_scaffold[2] >= 0)):
            if ((current_scaffold.initial_strand >= 0) or (current_scaffold.final_strand >= 0)):
                uid = uid + 1
                helix_topology[uid,0] = uid+1
                # set helix ID, lattice ID, and strand type: 0 for scaffold, 1 for staple
                helix_topology[uid,1] = num
                helix_topology[uid,2] = (i-1) + 1  # use 1-based IDs
                helix_topology[uid,3] = StrandType.SCAFFOLD
                uid_curr_scaf = helix_topology[uid,1:4]

                #  5'-neighbor
                helix_topology[uid,4] = current_scaffold.initial_strand
                helix_topology[uid,5] = current_scaffold.initial_base
                helix_topology[uid,6] = StrandType.SCAFFOLD

                # 3'-neighbor
                helix_topology[uid,7] = current_scaffold.final_strand
                helix_topology[uid,8] = current_scaffold.final_base  
                helix_topology[uid,9] = StrandType.SCAFFOLD

                # Watson-Crick neighbor
                #if ( (current_staple[0] >= 0) or (current_staple[2] >= 0)):
                if ( (current_staple.initial_strand >= 0) or (current_staple.final_strand >= 0)):
                    helix_topology[uid,10] = num
                    helix_topology[uid,11] = (i-1) + 1  # use 1-based IDs
                    helix_topology[uid,12] = 1
                else:
                    helix_topology[uid,10] = -1
                    helix_topology[uid,11] = -1
                    helix_topology[uid,12] = 1

                # coordinate
                helix_topology[uid,13] = scaffold_coords[i,0]
                helix_topology[uid,14] = scaffold_coords[i,1]
                helix_topology[uid,15] = scaffold_coords[i,2]

                # deletion
                helix_topology[uid,16] = deletions[i]

                # insertion
                helix_topology[uid,17] = insertions[i]
            #__if ((current_scaffold[0] >= 0) or (current_scaffold[2] >= 0))__

            # if the base exists in the staple strand
            if ((current_staple.initial_strand >= 0) or (current_staple.final_strand >= 0)):
                uid = uid + 1
                helix_topology[uid,0] = uid+1
                # set helix ID, lattice ID, and strand type: 0 for scaffold, 1 for staple
                helix_topology[uid,1] = num
                helix_topology[uid,2] = (i-1) + 1  # use 1-based IDs 
                helix_topology[uid,3] = StrandType.STAPLE
                uid_curr_stap = helix_topology[uid,1:4]

                #  5'-neighbor
                helix_topology[uid,4] = current_staple.initial_strand
                helix_topology[uid,5] = current_staple.initial_base  
                helix_topology[uid,6] = StrandType.STAPLE

                # 3'-neighbor
                helix_topology[uid,7] = current_staple.final_strand
                helix_topology[uid,8] = current_staple.final_base  
                helix_topology[uid,9] = 1 

                # Watson-Crick neighbor
                if ( (current_scaffold.initial_strand >= 0) or (current_scaffold.final_strand >= 0)):
                    helix_topology[uid,10] = num
                    helix_topology[uid,11] = (i-1) + 1
                    helix_topology[uid,12] = 0
                else:
                    helix_topology[uid,10] = -1
                    helix_topology[uid,11] = -1
                    helix_topology[uid,12] = 0

                # coordinate
                helix_topology[uid,13] = staple_coords[i,0]
                helix_topology[uid,14] = staple_coords[i,1]
                helix_topology[uid,15] = staple_coords[i,2]

                # deletion
                helix_topology[uid,16] = deletions[i]

                # insertion
                helix_topology[uid,17] = insertions[i]
            #__if ((current_staple[0] >= 0) or (current_staple[2] >= 0))__

            # check if a basepair exists
            if ( ((current_scaffold.initial_strand >= 0) or (current_scaffold.final_strand >= 0)) and 
                 ((current_staple.initial_strand >= 0) or (current_staple.final_strand >=0 ))):
                dnode = np.concatenate((dnode, dnode_0[[i],:]), axis=0)
                triad = np.concatenate((triad, triad_0[:,:,[i]]), axis=2)
                ida = np.concatenate((uid_curr_scaf, uid_curr_stap))
                id_nt_0 = np.concatenate((id_nt_0, [ida]), axis=0)
        #__for i in xrange(0,num_lattice)__

        # remove rows at the end that don't contain any information.
        helix_topology = np.delete(helix_topology, np.s_[uid+1:2*num_bases:1],0)
        return helix_topology, dnode, triad, id_nt_0, dnode_0

    def _generate_coordinates(self, lattice_type, row, col, strand_num, num_lattice):
        """Generate the coordinates for base pair positions along a helix and atom
           positions along the dna helix.
        """
        #self._logger.debug("-------------------- generate_coordinates --------------------")
        r_strand = 1.25      # half the distance between the axes of two neighboring DNA helices
        r_helix = 1.0        # radius of DNA helices (nm)
        dist_bp = 0.34       # rise between two neighboring base-pairs (nm)
        ang_bp = 360/10.5    # twisting angle between two neighboring base-pairs (degree)
        ang_minor = 120      # angle of the minor groove (degree)

        # positions of the scaffold nucleotide and staple nucleotide
        # in the local reference frame.
        scaf_local = r_helix * np.array([cos(180-ang_minor/2), sin(180-ang_minor/2), 0.0]).transpose();
        stap_local = r_helix * np.array([cos(180+ang_minor/2), sin(180+ang_minor/2), 0.0]).transpose();

        if (lattice_type == CadnanoLatticeType.honeycomb):
            xpos =  sqrt(3.0) * col * r_strand
            zpos = -3.0 * row * r_strand;

            if ( ((row % 2 == 0) and (col % 2 == 0)) or ((row % 2 != 0) and (col % 2 != 0))):
                zpos = zpos + r_strand

            if ( strand_num % 2 == 0 ):
                init_ang = -30 + ang_bp / 2
            else:
                init_ang = 150 + ang_bp / 2

        elif (lattice_type == CadnanoLatticeType.square): 
            xpos =  2.0 * col * r_strand
            zpos = -2.0 * row * r_strand
    
            if ( strand_num % 2 == 0):
                init_ang = 180 + ang_bp / 2
            else:
                init_ang = 0 + ang_bp / 2

        # Compute coordinate and twisting angle for each base pair. 
        init_coord_ang_strand = np.array([xpos, 0.0, zpos, init_ang]);
        scaffold_coords = np.zeros((num_lattice,3), dtype=float)
        staple_coords = np.zeros((num_lattice,3), dtype=float)
        dnode = np.zeros((num_lattice,3), dtype=float);
        triad = np.zeros((3,3,num_lattice), dtype=float);

        if (strand_num % 2 == 0):
            e3 = np.array([0, 1, 0],dtype=float)
        else:
            e3 = np.array([0, -1, 0],dtype=float)

        for i in xrange(0,num_lattice):
            ref_coord = init_coord_ang_strand + np.array([0, dist_bp*i, 0, ang_bp*i]);

            # Basepair position
            dnode[i,:] = ref_coord[0:3];

            # Basepair orientation
            e2 = np.array([cos(-_deg2rad(ref_coord[3])), 0, sin(-_deg2rad(ref_coord[3]))])
            e1 = np.cross(e2, e3);
            triad[:,:,i] = np.array([e1, e2, e3]).transpose();
            #print(">>> triad=%s" % str(triad[:,:,i]))

            # Scaffold & staple nucleotide positions
            scaffold_coords[i,:] = dnode[i,:] + np.dot(triad[:,:,i],scaf_local)
            staple_coords[i,:]   = dnode[i,:] + np.dot(triad[:,:,i],stap_local)

            #print(">>> ref_coord= %f %f %f %f" % (ref_coord[0], ref_coord[1], ref_coord[2], ref_coord[3]) )
            #logger.debug("scaffold_coords[%d,:]= %f %f %f" % (i+1, scaffold_coords[i,0],scaffold_coords[i,1],scaffold_coords[i,2]) )
            #print(">>> staple_coords[%d,:]= %f %f %f" % (i+1, staple_coords[i,0],staple_coords[i,1],staple_coords[i,2]) )
            #print(" ")
        #__for i in xrange(0,num_lattice)__

        return scaffold_coords, staple_coords, dnode, triad

    def _calculate_staple_ends(self, dna_topology):
        dna_topology,strands = self._build_strands(dna_topology)
        num_strands = len(strands)
        staple_ends = np.empty(shape=(0,5),dtype=float)

        for i in xrange(0,num_strands):
            strand = strands[i]
            is_scaf_strand = dna_topology[strand.tour[0]-1].is_scaf
            for j in xrange(0,len(strand.tour)):
                is_scaf_base = dna_topology[strand.tour[j]-1].is_scaf
            if (not is_scaf_strand):
                h0 = dna_topology[strand.tour[0]-1].h
                p0 = dna_topology[strand.tour[0]-1].p
                h1 = dna_topology[strand.tour[-1]-1].h
                p1 = dna_topology[strand.tour[-1]-1].p
                staple_ends = np.concatenate((staple_ends, [[i+1, h0, p0, h1, p1]]), axis=0)
        #__for i in xrange(0,num_strands)__
        return staple_ends

    def _build_strands(self,dna_topology):
        """ Build strands. 
        """
        num_bases = len(dna_topology)
        #print "[build_strand] num bases=%d" % num_bases 
        strands = []
        n_strand = 0
        is_visited = [False]*num_bases
        
        while (True):
            try:
                base_index = is_visited.index(False)
                curr_base = dna_topology[base_index]
                #print "[build_strand] base_index=%d" % base_index 
                #print "[build_strand] curr_base.id=%d" % curr_base.id 
            except:
                break

            init_baseID = curr_base.id
            strand = DnaStrand(n_strand)

            # Find the first base in the current strand
            while ((curr_base.up >= 0) and (curr_base.up != init_baseID)):
                curr_base = dna_topology[curr_base.up-1]
                if (is_visited[curr_base.id-1]):
                    sys.stderr.write('[build_strand] **** ERROR: Reached a visited base.');
                    return None,None
            #_while 

            strand = DnaStrand(n_strand)
            strands.append(strand)

            if (curr_base.up < 0):                     # currBase is at the 5'-end of the strand
                strand.is_circular = False
            elif (curr_base.up == init_baseID):        # currBase goes back to the starting point
                strand.is_circular = True
                curr_base = dna_topology[init_baseID-1]
            else:
                sys.stderr.write('[build_strand] **** ERROR: Exception.')
                return None,None

            # Walk through the current strand
            n_residue = 1
            strand.tour.append(curr_base.id)
            curr_base.strand = n_strand
            curr_base.residue = n_residue
            is_visited[curr_base.id-1] = True

            # Each loop adds a new base
            while ( (not strand.is_circular and (curr_base.down >= 0)) or 
                    (strand.is_circular and (curr_base.down != init_baseID)) ):
                curr_base = dna_topology[curr_base.down-1]

                if (is_visited[curr_base.id-1]):
                    sys.stderr.write('[build_strand] **** ERROR: Reached a visited base.\n')
                    return None,None

                if (n_residue == 1):
                    strand.is_scaffold = curr_base.is_scaf

                n_residue = n_residue + 1
                strand.tour.append(curr_base.id)
                curr_base.strand = n_strand
                curr_base.residue = n_residue
                is_visited[curr_base.id-1] = True
            #__while((not strand.is_circular 

            n_strand += 1
        #__while (True):
        return dna_topology, strands

    def set_sequence_from_name(self, modified_structure, seq_name):
        """ Set the sequence information for the staple and scaffold strands using a known
            origami vector sequence name.

            Cadnano seems to start assigning the sequence at virtual helix 0. If ordered_traverse=True
            then use the cadnano method.

            If the structure has not been modified with insertions and deletions then we will need to 
            use its insertions and deletions information (at the base level) to selectively set its 
            sequence.

            Arguments:
                modified_structure (bool): If true then the structure has been modified for insertions and deletions. 
                seq_name (string): The name of the sequence as defined in the dna_sequence_data dictionary. 
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.setLevel(logging.INFO)

        sequence = dna_sequence_data.get(seq_name,None)
        seq_index = 0
        sequence_length = len(sequence)
        self._logger.debug("-------------------- set_sequence_from_name --------------------")
        self._logger.debug("sequence name: %s" % seq_name)
        self._logger.debug("sequence length: %d" % len(sequence))
        dna_structure = self.dna_structure 
        base_connectivity = dna_structure.base_connectivity
        strands = dna_structure.strands 
        ordered_traverse = False
        ordered_traverse = True

        for strand in strands:
            if (not strand.is_scaffold):
                continue 

            if (ordered_traverse): 
                self._logger.debug("---------- traverse scaffold strand --------------------")
                tour = strand.tour
                id = tour[0]
                base = base_connectivity[id-1]
                min_vh = base.h
                min_p = base.p
                start_index = 0

                for j in xrange(0,len(tour)):
                    id = tour[j]
                    base = base_connectivity[id-1]
                    if (base.h < min_vh):
                        min_vh = base.h
                        min_p = base.p
                        start_index = j
                    if ((base.h == min_vh) and (base.p < min_p)):
                        min_p = base.p
                        start_index = j
                self._logger.debug("start_index %d  min_vh %d  min_p %d ", start_index, min_vh,min_p)
                start_id = tour[start_index]
                start_base = base_connectivity[start_id-1]
                base = start_base 
                even_vh = min_vh % 2

                for j in xrange(0,len(tour)):
                    letter = sequence[seq_index]
                    up = base.up
                    down = base.down
                    across = base.across

                    if (not modified_structure):
                        if (across >= 0):
                            if (base.skip != 0):
                                self._logger.debug("**** deleted base: id %d " % base.id)
                                letter = 'N'
                            elif (base.loop != 0):
                                self._logger.debug("**** inserted base: id %d " % base.id)
                                strand.insert_seq.append(seq.letters[seq_index+1])
                                seq_index += 2
                            else:
                                seq_index += 1
                    else:
                        seq_index += 1

                    if (seq_index == sequence_length): 
                        seq_index = 0

                    base.seq = letter 
                    self._logger.debug("base id %d  vh %d  pos %d  up %d  down %d  across %d  seq %s", 
                        base.id, base.h, base.p, up, down, across, base.seq)

                    if (across >= 0):
                        across_base = base_connectivity[base.across-1]
                        across_base.seq = self._wspair(letter)

                    if (even_vh):
                        if (up != -1):
                            base_id = up
                            base = base_connectivity[base_id-1]
                        elif (across != -1):
                            base_id = across
                            base = base_connectivity[base_id-1]
                            even_vh = base.h % 2
                    else:
                        if (down != -1):
                            base_id = down
                            base = base_connectivity[base_id-1]
                        elif (across != -1):
                            base_id = across
                            base = base_connectivity[base_id-1]
                            even_vh = base.h % 2

                 #__for j in xrange(0,len(tour)):

            else: 
                for i in xrange(0,len(strand.tour)):
                    letter = sequence[seq_index]
                    base_index = int(strand.tour[i])-1
                    base = base_connectivity[base_index]
                    across = base.across

                    if (not modified_structure):
                        if (base.skip != 0):
                            self._logger.debug("**** deleted base: id %d  vh %d  pos %d", base.id, base.h, base.p)
                            letter = 'N'
                        elif (base.loop != 0):
                            self._logger.debug("**** inserted base: id %d  vh %d  pos %d " % base.id, base.h, base.p)
                            strand.insert_seq.append(seq.letters[seq_index+1])
                            seq_index += 2
                        else:
                            seq_index += 1
                    else:
                        seq_index += 1

                    base.seq = letter 

                    if (seq_index == sequence_length): 
                        seq_index = 0

                    if (base.across >= 0):
                        across_base = base_connectivity[base.across-1]
                        across_base.seq = self._wspair(letter)

                #__for i in xrange(0,len(strand.tour))
        #__for strand in strands

        print_strands = True
        print_strands = False
        if (print_strands):
            self._logger.debug("---------- strands sequences ----------")
            self._logger.debug(">>> number of strands: %d " % len(strands))
            for strand in strands:
                tour = strand.tour 
                self._logger.debug(">>> strand: %d  scaf: %d length: %d" % (strand.id,strand.is_scaffold,len(tour)))
                self._logger.debugprint("    seq:")
                for j in xrange(0,len(tour)):
                    id = tour[j]
                    base = base_connectivity[id-1]
                    self._logger.debug("    vhelix: %d  pos: %d  seq: %s" % (int(base.h), int(base.p), base.seq))
            #__for i in xrange(0,len(strands))

    #__def set_sequence_from_name

    def set_sequence(self, modified_structure, sequence):
        """ Set the sequence information for the staple and scaffold strands.

            The caDNAno csv file contains sequence information reflecting the insertions and deletions of a
            design. If a structure has had bases inserted and deleted then its sequence can be set directly 
            from the sequences in the csv file. If the structure has not been modified with insertions and
            deletions then we will need to use its insertions and deletions information (at the base level)
            to set selectively set its sequence from the sequences in the csv file.   

            Arguments:
               modified_structure (bool): If True then the structure has been modified with deletions and insertions. 
               sequence (DnaSequence): A list of DnaSequence objects representing the sequences for staple or scaffold strands.
        """

        #self._logger.setLevel(logging.DEBUG)
        self._logger.setLevel(logging.INFO)

        dna_structure = self.dna_structure 
        strands = dna_structure.strands 
        base_connectivity = dna_structure.base_connectivity
        staple_ends = dna_structure.staple_ends 
        start = np.array([0,0], dtype=float)

        for i in xrange(0,len(sequence)):
            seq = sequence[i]
            start[0] = sequence[i].start[0]
            start[1] = sequence[i].start[1]
            row = _find_row(start, staple_ends[:,1:3])[0]
            istrand = int(staple_ends[row,0])
            tour = strands[istrand-1].tour

            if (modified_structure):
                for j in xrange(0,len(strands[istrand-1].tour)):
                    k = int(strands[istrand-1].tour[j])-1
                    base_connectivity[k].seq = seq.letters[j]
                    if (self.base_connectivity[k].across >= 0):
                        k_across = self.base_connectivity[k].across
                        self.base_connectivity[k_across-1].seq = self._wspair(seq.letters[j])
                #__for j

            else:
                seq_index = 0
                for j in xrange(0,len(tour)):
                    letter = seq.letters[seq_index]
                    base_index = int(tour[j])-1
                    base = base_connectivity[base_index]
                    if (base.skip != 0):
                        letter = 'N'
                    elif (base.loop != 0):
                        strand.insert_seq.append(seq.letters[seq_index+1])
                        seq_index += 2
                    else:
                        seq_index += 1
                    base.seq = letter
                    if (base.across >= 0):
                        across_base = base_connectivity[base.across-1]
                        across_base.seq = self._wspair(letter)
                #__for j

        #__for i

        print_strands = True
        print_strands = False
        if (print_strands):
            self._logger.debug("---------- strands sequences ----------")
            self._logger.debug(">>> number of strands: %d " % len(strands))
            for strand in strands:
                tour = strand.tour
                self._logger.debug(">>> strand: %d scaf: %d len: %d" % (strand.id,strand.is_scaffold,len(tour)))
                self._logger.debug("    seq:")
                for j in xrange(0,len(tour)):
                    id = tour[j]
                    base = base_connectivity[id-1]
                    self._logger.debug("    vhelix: %d  pos: %d  seq: %s" % (int(base.h), int(base.p), base.seq))
            #__for i in xrange(0,len(strands))
    #__def set_sequence

    def _wspair(self, x):
        """ Match a base with its complementary base. 
        """
        x = x.upper()

        if (x == 'A'):
            y = 'T'
        elif (x == 'G'):
            y = 'C'
        elif (x == 'C'):
            y = 'G'
        elif (x == 'T'):
            y = 'A'
        elif (x == 'N'):
            y = 'N'
        else:
            self._logger.error('Illegal base.')
        return y

#__class CadnanoTopology(object)

class _Timer(object):
    """ The Timer class is used to calculare elapsed time between calls to the start and finish methods."""
    def __init__(self):
        self.start_time = 0.0
        self.end_time = 0
        self.secs = 0
        self.msecs = 0

    def start(self):
        self.start_time = time.time()

    def finish(self):
        self.end_time = time.time()
        self.secs = self.end_time - self.start_time
        self.msecs = self.secs * 1000  # millisecs
        return self.secs

def _deg2rad(deg):
    """Convert degrees into radians. """
    rad = (pi/180)* deg;
    return rad

def _find_row(neigh, curr_bases):
    tmp = curr_bases - neigh
    s = np.sum(np.abs(tmp),1)
    ind = np.where(s==0)[0]
    return ind



def _vrrotmat2vec(R):
    """ Extract the equivalent rotation about an axis from a rotation matrix. """
    m00 = R[0,0]
    m01 = R[0,1]
    m02 = R[0,2]
    m10 = R[1,0]
    m11 = R[1,1]
    m12 = R[1,2]
    m20 = R[2,0]
    m21 = R[2,1]
    m22 = R[2,2]
    angle = acos(( m00 + m11 + m22 - 1)/2.0)
    x = (m21 - m12) / sqrt( pow(m21-m12,2) + pow(m02-m20,2) + pow(m10-m01,2) )
    y = (m02 - m20) / sqrt( pow(m21-m12,2) + pow(m02-m20,2) + pow(m10-m01,2) )
    z = (m10 - m01) / sqrt( pow(m21-m12,2) + pow(m02-m20,2) + pow(m10-m01,2) )
    return np.array([x,y,z],dtype=float),angle


def _vrrotvec2mat(axis, theta):
    """ Create a rotation matrix to rotate theta degrees about the axis defined by vec. """
    s = np.sin(theta)
    c = np.cos(theta)
    t = 1 - c
    #print("[vrrotvec2mat] theta=%s" % str(theta))

    # normalize the vector
    x,y,z = axis / np.linalg.norm(axis)
    #print("[vrrotvec2mat] x=%s" % str(x))
    #print("[vrrotvec2mat] y=%s" % str(y))
    #print("[vrrotvec2mat] z=%s" % str(z))

    return np.array([ [t*x*x + c,   t*x*y - s*z,  t*x*z + s*y],
                      [t*x*y + s*z, t*y*y + c,    t*y*z - s*x],
                      [t*x*z - s*y, t*y*z + s*x,  t*z*z + c  ]])


def _bp_interp(dnode_1, triad_1, dnode_2, triad_2, n):
    """ Interpolate the position (dnode) and orientation (triad) between two base pairs.
        Solve for rotation matrix R defined as the solution of: 
            R * triad_1 = triad_2 (multiply both sides by transpose(triad_1).
    """
    #print("[bp_interp] n=%s " % str(n))
    dnode_interp = np.zeros((n,3),dtype=float)
    triad_interp = np.zeros((3,3,n),dtype=float)

    # solve for rotation matrix R
    R = np.dot(triad_2,triad_1.T)
    #print("[bp_interp] R=%s " % str(R))
    #print("[bp_interp] R.shape=%s " % str(R.shape))

    a,theta = _vrrotmat2vec(R)
    #print("[bp_interp] a=%s " % str(a))
    #print("[bp_interp] theta=%s " % str(theta))

    # Calculate for dNode_interp and triad_interp
    for i in xrange(0,n):
        dnode_interp[i,:] = (dnode_1*(n+1-(i+1)) + dnode_2*(i+1)) / (n+1)
        angle = theta*(i+1)/(n+1)
        #print("[bp_interp] angle=%s " % str(angle))
        rot_mat = _vrrotvec2mat(a, angle)
        #print("[bp_interp] rot_mat=%s\n" % str(rot_mat))
        #print("[bp_interp] rot_mat.shape=%s\n" % str(rot_mat.shape))
        triad_interp[:,:,i] = np.dot(rot_mat,triad_1)

    return dnode_interp, triad_interp


def main():
    """ Create a topology table from a caDNAno JSON design file."""
    json_file_name = sys.argv[1]
    cadnano_reader = CadnanoReader()
    #cadnano_reader.set_logging_level(logging.DEBUG)
    cadnano_model = cadnano_reader.read_json(json_file_name)
    convert_design = CadnanoConvertDesign()
    dna_structure = convert_design.create_structure(cadnano_model)

    if (len(sys.argv) == 3):
       csv_file_name = sys.argv[2]
       seq = cadnano_reader.read_csv(csv_file_name)
       structure.set_sequence(seq)

if __name__ == '__main__':
    main()

