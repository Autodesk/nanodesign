#!/bin/sh -e

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

# This file is a shortcut to run the tests and produce the JUnitXML output to
# stdout. This is needed for Jenkins CI automation of the tests.

# Don't assume this script is run in same directory as where this file and the test file reside. 
# Optionally set OUTPUT_DIR (with trailing slash) to location of output file
pathPrefix=$(dirname $0)

echo "python -m pytest --tb=short $pathPrefix/tests_basic.py --junitxml=${OUTPUT_DIR}output.xml"
python -m pytest --tb=short $pathPrefix/tests_basic.py --junitxml=${OUTPUT_DIR}output.xml
