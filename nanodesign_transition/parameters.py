#!/usr/bin/env python
"""
This module is used to store information for DNA physical parameters.
"""

class DnaParameters:
    """ This class stores information for DNA parameters.

        Attributes:
            helix_radius: The radius of DNA helices (nm).
            strand_radius: The half the distance between the axes of two neighboring DNA helices.
            base_pair_rise: The rise between two neighboring base-pairs (nm).
            base_pair_twist_angle: The twisting angle between two neighboring base-pairs (degree).
            minor_groove_angle: The angle of the minor groove (degree).
    """
    helix_radius = 1.0
    strand_radius = 1.25
    base_pair_rise = 0.34
    base_pair_twist_angle = 360 / 10.5
    minor_groove_angle = 120.0

