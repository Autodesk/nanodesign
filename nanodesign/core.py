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
