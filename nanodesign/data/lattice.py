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
This module contains classes used to define the lattice-base geometric framework of a DNA structure.

"""
from abc import ABCMeta, abstractmethod
import numpy as np
from math import sqrt

from .parameters import DnaParameters
from ..converters.cadnano.common import CadnanoLatticeType

class Lattice(object):
    """ This is the lattice base class.
    """
    __metaclass__ = ABCMeta

    # Set the radius using a physical value or one used for visualization.
    def __init__(self, radius):
        self.radius = radius

    @abstractmethod
    def get_neighbors(row,col):
        """ Get the neighboring lattice coordinates for a given lattice coordinate. """
        raise NotImplementedError

    @abstractmethod
    def get_neighbor_direction(row, col, nrow, ncol):
        """ Get the direction for a neighboring lattice coordinate. """
        raise NotImplementedError

    @abstractmethod
    def get_neighbor_index(row, col, nrow, ncol):
        """ Get the index for a neighboring lattice coordinate. """
        raise NotImplementedError

    @staticmethod
    def create_lattice(lattice_type, helix_distance):
        """ Create a lattice object for a given type. 

            Arguments:
                lattice_type (Int): The lattice type (e.g. honeycomb) taken from CadnanoLatticeType.
                helix_distance (Float): The distance between adjacent DNA helices.
        """
        radius = helix_distance / 2.0
        if (lattice_type == CadnanoLatticeType.honeycomb):
            return HoneycombLattice(radius)
        elif (lattice_type == CadnanoLatticeType.square):
            return SquareLattice(radius)
        else:
            return None

    @staticmethod
    def even_parity_coordinate(row, column):
        """ Determine if the given lattice coordinate is even parity. """ 
        return (row % 2) == (column % 2)

    @staticmethod
    def odd_parity_coordinate(row, column):
        """ Determine if the given lattice coordinate is odd parity. """ 
        return (row % 2) ^ (column % 2)


class SquareLattice(Lattice):
    """ This class defines the data and methods for a square lattice.
        
        Attributes: 
            index (List[Tuple()]): The list used to map lattice (row,col) coordinates to a neighbor index (0-3). 
            number_of_neighbors (int): The number of lattice cell neighbors.
            scaffold_low (List[List[Int]]): The list scaffold of crossover points.
            scaffold_high (List[List[Int]]): The list scaffold of crossover points.
            staple_low (List[List[Int]]): The list staple of crossover points.
            staple_high (List[List[Int]]): The list staple of crossover points.
            step (int): The number of bases between crossover positions. 

        Some of the data in this class is used to compute crossover positions between neighboring helices on a square lattice. 
        It is based on the data from the caDNAno squarepart.py file. 

        The index_map list is used to map a lattice (row,col) coordinate to a neighbor (nrow,ncol) index. The index 
        identifies the neighbor direction in a standard way.

    """
    def __init(self,radius):
       Lattice.__init__(self, radius)
    number_of_neighbors = 4
    step = 32
    scaffold_low  = [ [4, 26, 15], [18, 28, 7], [10, 20, 31], [2, 12, 23] ]
    scaffold_high = [ [5, 27, 16], [19, 29, 8], [11, 21, 0], [3, 13, 24]  ]
    staple_low    = [ [31], [23], [15], [7] ]
    staple_high   = [  [0], [24], [16], [8] ]
    index_map = [ (0,1), (-1,0), (0,-1), (1,0) ]

    @classmethod
    def get_neighbors(cls,row,col):
        """ Get the neighboring lattice coordinates for a given lattice coordinate.

            Arguments:
                row (int): The row lattice coordinate. 
                col (int): The column lattice coordinate. 

            Returns:
                neighbors (list[(int,int)): The list of neighboring lattice coordinates. 
        """
        neighbors = []

        if Lattice.even_parity_coordinate(row, col):
            neighbors.append((row,   col+1))
            neighbors.append((row+1, col  ))
            neighbors.append((row,   col-1))
            neighbors.append((row-1, col  ))
        else:
            neighbors.append((row,   col-1))
            neighbors.append((row-1, col  ))
            neighbors.append((row,   col+1))
            neighbors.append((row+1, col  ))

        return neighbors

    @classmethod
    def get_neighbor_direction(cls, row, col, nrow, ncol):
        """ Get the direction for a neighboring lattice coordinate. """
        dx = nrow - row 
        dz = col - ncol
        return [dx, dz]

    @classmethod
    def get_neighbor_index(cls, row, col, nrow, ncol):
        """ Get the index for a neighboring lattice coordinate. 
            Arguments:
                row (int): The row lattice coordinate.
                col (int): The column lattice coordinate.
                nrow (int): The neighbor row lattice coordinate.
                ncol (int): The neighbor column lattice coordinate.

            Returns:
               nindex (Int): The neighbor index or None if the lattice coordinates are not a neighbor. 
        """
        nval = (nrow-row, ncol-col)
        if nval in cls.index_map:
            nindex = cls.index_map.index(nval)
        else:
            nindex = None
        return nindex 

class HoneycombLattice(Lattice):
    """ This class defines the data and methods for a honeycomb lattice.

        Attributes: 
            index_map_odd (List[Tuple()]): The list used to map lattice (row,col) coordinates into a neighbor index (0-2)
               for odd-numbered helices. 
            index_map_even (List[Tuple()]): The list used to map lattice (row,col) coordinates into a neighbor index (0-2) 
               for even-numbered helices. 
            number_of_neighbors (int): The number of lattice cell neighbors.
            scaffold_low (List[List[Int]]): The list scaffold of crossover points.
            scaffold_high (List[List[Int]]): The list scaffold of crossover points.
            staple_low (List[List[Int]]): The list staple of crossover points.
            staple_high (List[List[Int]]): The list staple of crossover points.
            step (int): The number of bases between cross-over positions (step size). 

        Some of the data in this class is used to compute crossover positions between neighboring helices 
        on a honeycomb lattice. It is based on the data from the caDNAno honeycombpart.py file. 

        The index_map list is used to map a lattice (row,col) coordinate to a neighbor (nrow,ncol) index. The index 
        identifies the neighbor direction in a standard way.
    """
    def __init(self, radius):
       Lattice.__init__(self, radius)
    number_of_neighbors = 3
    step = 21
    scaffold_low  = [ [1, 11], [8, 18], [4, 15] ]
    scaffold_high = [ [2, 12], [9, 19], [5, 16] ]
    staple_low    = [ [6], [13], [20] ]
    staple_high   = [ [7], [14], [0]  ]
    index_map_odd  = [ ( 1,0), (0,1), (0,-1) ]
    index_map_even = [ (-1,0), (0,1), (0,-1) ]

    @classmethod
    def get_neighbors(cls,row,col):
        """ Get the neighboring lattice coordinates for a given lattice coordinate.

            Arguments:
                row (int): The row lattice coordinate. 
                col (int): The column lattice coordinate. 

            Returns:
                neighbors (list[(int,int)): The list of neighboring lattice coordinates. 
        """
        neighbors = []

        if Lattice.even_parity_coordinate(row, col):
            neighbors.append((row,  col+1))
            neighbors.append((row-1,col))
            neighbors.append((row,  col-1))
        else:
            neighbors.append((row,   col-1))
            neighbors.append((row+1, col  ))
            neighbors.append((row,   col+1))

        return neighbors

    @classmethod
    def get_neighbor_direction(cls, row, col, nrow, ncol):
        """ Get the direction for a neighboring lattice coordinate. 
        """
        # TODO(DaveP) We can do this using a lookup table. 
        radius = Lattice.radius
        x = sqrt(3.0)* radius * col
        z =     -3.0 * radius * row
        if Lattice.odd_parity_coordinate(row, col):
            z += radius

        nx = sqrt(3.0) * radius * ncol
        nz =     -3.0  * radius * nrow
        if Lattice.odd_parity_coordinate(nrow, ncol):
            nz += radius

        dx = nx - x
        dz = nz - z
        mag = sqrt(dx*dx + dz*dz)

        if mag != 0.0:
            dx = dx / mag
            dz = dz / mag
        else:
            dx = 0.0
            dz = 0.0

        return [dx,dz]

    @classmethod
    def get_neighbor_index(cls, row, col, nrow, ncol):
        """ Get the index for a neighboring lattice coordinate.

            Arguments:
                row (int): The row lattice coordinate.
                col (int): The column lattice coordinate.
                nrow (int): The neighbor row lattice coordinate.
                ncol (int): The neighbor column lattice coordinate.

            Returns:
               nindex (Int): The neighbor index or None if the lattice coordinates are not a neighbor.   
        """
        nval = (nrow-row, ncol-col)
        # Get the odd or even map depending on the (row,col).
        if cls.even_parity_coordinate(row, col):
           index_map = cls.index_map_even
        else:
           index_map = cls.index_map_odd  
        # Map the neighbor coordinate to an index.
        if nval in index_map:
            nindex = index_map.index(nval)
        else:
            nindex = None
        return nindex 

