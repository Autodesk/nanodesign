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
The (dummy) InternalData class. Probably we can remove this and replace with a better name, etc, but this gives a basic skeleton for what it might look like.
"""

__all__ = ['InternalData']

# please use the new style classes!
class InternalData(object):
    """ The internal data type used by most of the algorithms. """
    pass

class DummyHelper(object):
    """ This is used internally, but shouldn't be exposed as part of the module. Note how __all__ doesn't include this name. """
    pass


