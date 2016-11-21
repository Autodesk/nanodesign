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
This module is used to store information for DNA physical parameters.
"""

class DnaParameters(object):
    """ This class stores information for DNA parameters.

        Attributes:
            helix_distance: The distance between the axes of two neighboring DNA helices.
            helix_radius: The radius of DNA helices (nm).
            base_pair_rise: The rise between two neighboring base-pairs (nm).
            base_pair_twist_angle: The twisting angle between two neighboring base-pairs (degree).
            minor_groove_angle: The angle of the minor groove (degree).
    """
    def __init__(self):
        self.helix_distance = 2.3
        self.helix_radius = 1.0
        self.base_pair_rise = 0.34
        self.base_pair_twist_angle = 360 / 10.5
        self.minor_groove_angle = 120.0

class DnaPolarity:
    """ This class defines constants used to describe the polarity of a DNA strand. 

        Attributes:
            FIVE_PRIME (string): The direction of a DNA strand is from 5'->3'. When indexing the bases in a helix this means 
                indexes increase in the 5'->3' direction.
            THREE_PRIME (string): The direction of a DNA strand is from 3'->5'. When indexing the bases in a helix this means 
                indexes increase in the 3'->5' direction.
    """
    FIVE_PRIME = "5'"
    THREE_PRIME = "3'"

class DnaBaseNames:
    A = "A"
    C = "C"
    G = "G"
    T = "T"

