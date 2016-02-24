##############################
# COPYRIGHT notice goes here #
##############################


# This file is the basic package setup. It's currently fleshed out with a very minimal package list. If this goes into the PyPI repository, we need to add appropriate metadata, etc.

from distutils.core import setup

if __name__ == '__main__':
    setup(name="nanodesign", version="1.0",
          packages=['nanodesign','nanodesign.converters','nanodesign.algorithms','nanodesign.data'],
      )
