#!/bin/bash

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

python -m pytest --tb=short tests_basic.py --junitxml=output.xml >> /dev/null
export code=$?
cat output.xml
exit $code

