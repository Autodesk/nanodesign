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
   core routines that are part of the nanodesign base package

   Currently some sample code.
"""

__all__ = ['load']


def load(filename, type='auto'):
    """Load an arbitrary file which is assumed to be some nanostructure design file
    type, such as Cadnano, Cando, Tiamat, etc.
    
    Parameters
    ----------
    filename : str
        The filename to load; can include path specifiers, just like would be
        expected for Python's open command.
    type : {'auto' (default), 'cadnano', 'cando', 'tiamat'} optional
        Specifies the type of the file being loaded. If set to 'auto' it will
        try to figure it out and return an exception if it can't. Otherwise
        assumes it is of the type given and uses the appropriate converter.

    Returns
    -------
    filedata : data.InternalData
    
    Raises
    ------
    Some exceptions based on standard python open, plus specific ones for when type could not be detected. (more details should go here)
    """

    # Right now this is just going to return a string containing the first line of the file, rather than any real data.
    f = open(filename,'rt')
    filedata = f.readline()
    return filedata
