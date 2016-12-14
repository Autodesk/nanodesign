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
This module is used to store data for a DNA base. 

A list of DnaBase objects, stored in the DnaStructure object, is used to define the connectivity (topology) 
of a DNA structure: bases that are adjacent along a single DNA helix and which base they are paired with (if any). 
"""

class DnaBase(object):
    """ This class stores information for a DNA base. 

        Attributes:
            across (VisBase): The base's Watson-Crick neighbor.
            coordinates ((3x1 numpy float arrayList[Float]): The base helix axis coordinates.
            domain (int): The domain ID the base is in.
            down (VisBase): The base's 3' neighbor.
            h (int): The ID of the helix the base is in. 
            id (int): The base ID.
            is_scaf (bool): If True then this base is in a scaffold strand. 
            num_insertions (int): The number of insertions at this base.
            nt_coords ((3x1 numpy float arrayList[Float]): The base nucleotide coordinates.
            ref_frame ((3x1 numpy float arrayList[Float]): The base helix axis reference frame.
            p (int): The helix position of the base.
            seq (string): A one character string representing the base sequence nucleotide. 
            num_deletions (int): The number of deletions at this base.
            strand (int): The strand ID the base is in.
            up (VisBase): The base's 5' neighbor.

        The base coordinates and reference frame are references to elements of arrays stored in 
        the helix they are associated with.
    """

    def __init__( self, id, up=None, down=None, across=None, seq='N'):
        self.id = int(id)
        self.up = up
        self.down = down
        self.across = across
        self.seq = seq
        self.strand = None
        self.domain = None
        self.num_insertions = 0
        self.num_deletions = 0
        self.h = -1
        self.p = -1
        self.nt_coords = None
        self.coordinates = None
        self.ref_frame = None
        self.is_scaf = False

    def remove(self):
        """ Remove the base from the DNA structure.

            This function is used to remove bases flagged for deletion in the design file. 
            and is called when the design file is being processed. It will reset neighboring
            bases connectivity including the paired base if there is one. 
        """
        # Find the four neighboring bases.
        neighbor_up = self.up
        neighbor_down = self.down
        across = self.across
        if across != None:
            neighbor_across_up = across.up
            neighbor_across_down = across.down
        else:
            neighbor_across_up = None
            neighbor_across_down = None

        # Update base connectivity.
        if neighbor_up != None:
            neighbor_up.down = neighbor_down
        if neighbor_down != None:
            neighbor_down.up = neighbor_up

        # Update paired base connectivity.
        if neighbor_across_up != None:
            neighbor_across_up.down = neighbor_across_down
        if neighbor_across_down != None:
            neighbor_across_down.up = neighbor_across_up

    #__def remove

#__class DnaBase(object)


