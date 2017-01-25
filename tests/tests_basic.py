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

import pytest
import subprocess

import os.path
import hashlib

###################
# Setup path data #
###################

tests_path = os.path.dirname( os.path.abspath( __file__ ))
samples_path = os.path.join( tests_path, 'samples/')
base_path = os.path.abspath( os.path.join( tests_path, '../' ))
scripts_path = os.path.join( base_path, 'scripts/')

####################
# Helper Functions #
####################

def fast_hash_file( filename ):
    md5 = hashlib.md5()
    with open(filename, 'rb') as f:
        while True:
            data = f.read( 1048576 )  # 1 MB
            if not data:
                break
            md5.update( data )
    return md5.hexdigest()

# Note on master hashfile formats used here:
#
# We assume that the first entry, the filename, has no path information. It
# currently is assumed to reside in the tests/samples/ directory, but in the
# future we may make there be subdirectories there. If we do so, we will need to
# examine the functions and fixtures here that call os.path.basename, as well as
# the samples_path variable, and possibly more.

def load_master_hashfile(filename = "tests_master_hashfile.txt"):
    file_path = os.path.join( tests_path, filename )
    lines = []
    master = {}
    with open( file_path, 'rt') as f:
        lines = f.readlines()
        
    for line in lines:
        segments = line.split()
        name = segments[0]
        hashes = {}
        for segment in segments[1:]:
            key,value = segment.split(':')
            hashes[key] = value
        master[name] = hashes
    return master

def add_to_hashfile( hashfile, file_path, key, value ):
    filename = os.path.basename( file_path )    
    try:
        hashfile[filename][key] = value
    except KeyError:
        hashfile[filename] = {}
        hashfile[filename][key] = value
   

def save_master_hashfile():
    pass


#######################
# Setup Hashfile Data #
#######################

master_hashfile = load_master_hashfile()
current_tests_hashfile = {}

############
# Fixtures #
############


@pytest.fixture(scope="module",
                params=["fourhelix.json", "flat_sheet.json","beachball.json"])
def sample_file( request ):
    return os.path.join(samples_path, request.param)



@pytest.fixture()
def basic_converter_hash( sample_file ):
    filename = os.path.basename( sample_file ) 
    try:
        return master_hashfile[filename]['converter_basic']
    except KeyError:
        return None


#########
# Tests #
#########

def test_convert_basic( sample_file, basic_converter_hash ):
    converter_file = os.path.join( scripts_path, 'converter.py' )
    result = subprocess.call([converter_file,"--infile", sample_file , "--informat", "cadnano", "--inseqname", "M13mp18","--outfile", "my_sample_viewer.json", "--outformat", "viewer"] , stdout=None, stderr=None)
    # First way it could fail is if the call did not succeed, e.g. some error while executing.
    assert result == 0

    # Second way it could fail is if the md5 hash was different.
    # result = subprocess.check_output(['md5','-q','my_sample_viewer.json'])
    # result = result.rstrip('\n')
    # assert result == md5_hash

    result = fast_hash_file('my_sample_viewer.json')
    add_to_hashfile( current_tests_hashfile, sample_file, 'converter_basic', result )
    assert result == basic_converter_hash, "Hash value mismatch."


def test_convert_modify():
    filename = os.path.join( samples_path, 'flat_sheet.json' )
    converter_file = os.path.join( scripts_path, 'converter.py' )
    result = subprocess.call([converter_file,"--infile", filename , "--informat", "cadnano", "--inseqname", "M13mp18","--modify","true","--outfile", "my_sample_viewer.json", "--outformat", "viewer"] , stdout=None, stderr=None)
    # First way it could fail is if the call did not succeed, e.g. some error while executing.
    assert result == 0

    result = fast_hash_file('my_sample_viewer.json')
    add_to_hashfile( current_tests_hashfile, filename, 'converter_modify', result )
    assert result == master_hashfile['flat_sheet.json']['converter_modify'], "Hash value mismatch."



def test_show_hashes( capsys ):
    """This is a dummy test that should be run last. It will always fail, and will
show hash information in the stdout section of its failure report."""

    if capsys is not None:
        out,err = capsys.readouterr()

    import pprint
    pprint.pprint(current_tests_hashfile)

    # Change this line if you want to get a printout of all hashes processed.
    assert 0 == 0





    
