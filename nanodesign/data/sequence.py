#!/usr/bin/env python
"""
This module is used to store information for a DNA sequence. 
"""

class DnaSequence(object):
    def __init__(self, start, end, letters, length):
        self.start = list(start)
        self.end = list(end)
        self.letters = letters
        self.length = length
