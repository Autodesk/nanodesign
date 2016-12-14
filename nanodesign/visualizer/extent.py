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

""" This module is used to store coordinate extent information.  """

class VisExtent(object):
    """ This class stores extent information. """
    def __init__(self):
        self.is_set = False
        self.xmin = 0.0
        self.ymin = 0.0
        self.zmin = 0.0
        self.xmax = 1.0
        self.ymax = 1.0
        self.zmax = 1.0

    def set(self,xmin,xmax,ymin,ymax,zmin,zmax):
        """ Set the exetent with the given arguments. """
        self.xmin = xmin
        self.ymin = ymin
        self.zmin = zmin
        self.xmax = xmax
        self.ymax = ymax
        self.zmax = zmax
        self.is_set = True

    def reset(self):
        """ Reset the exetent to a unit cube. """
        self.is_set = False
        self.set(0.0,1.0,0.0,1.0,0.0,1.0)

    def update(self,x,y,z):
        """ Update the exetent with the given (x,y,z) coordinate. """

        if not self.is_set:
            self.set(x,x,y,y,z,z)
            return

        if (x < self.xmin):
            self.xmin = x
        elif (x > self.xmax):
            self.xmax = x

        if (y < self.ymin):
            self.ymin = y
        elif (y > self.ymax):
            self.ymax = y

        if (z < self.zmin):
            self.zmin = z
        elif (z > self.zmax):
            self.zmax = z

    def get_bounds(self):
        return self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax

    def get_center(self):
        return (self.xmin+self.xmax)/2.0, (self.ymin+self.ymax)/2.0, (self.zmin+self.zmax)/2.0

    def get_widths(self):
        return (self.xmax-self.xmin), (self.ymax-self.ymin), (self.zmax-self.zmin)
