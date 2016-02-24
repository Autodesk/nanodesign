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


