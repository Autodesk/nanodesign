Nanodesign
==========

Create, interact and modify nucleic acid based structures. The Nanodesign Python Package is aimed at providing a toolkit for working with structural DNA/RNA nanotechnology designs. This package is meant to support all types of interaction with these designs, including loading Cadnano files and modifying them, building new structures from scratch, or just converting formats so you can easily visualize or simulate your design.

## Getting Started

You should have [`git`](https://git-scm.com/downloads) installed as well as [Python 2.7](https://python.org). Note that this part of the guide is aimed towards Mac OS X or Linux installs. For a Windows install, a similar pattern should work using one of the command line shells.

On the command line, change to the directory you want Nanodesign installed under, and run:
```shell
git clone https://github.com/Autodesk/nanodesign
```

You will now have the repository cloned and checked out to the `master` branch in the subdirectory `nanodesign/`. See [CONTRIBUTING.md](CONTRIBUTING.md#branching-model) for details on other branches.

### Installation

We are currently working on getting the package set up for installation in your Python's site packages. Once this is fully tested, we'll also be adding the package to the PyPI package repository and it should no longer be necessary to acquire the git repository for regular usage. 

To try out the site package install, you can run the following from the repository directory:
```shell
python setup.py install
```

Currently, the install should be functional for almost all of the package, except for the [PDB/mmCIF](http://mmcif.wwpdb.org/) export routines. The example scripts will not be installed, but you can still access them from the `scripts/` subdirectory as mentioned in the examples, below.

### Example Usage: Package

To do some simple operations with the package, we have a few short scripts in the [`examples/`](examples/) subdirectory that show off the interface:

#### Strand Statistics

Source: [`examples/strand-statistics.py`](examples/strand-statistics.py)  
This script loads a Cadnano file and produces some statistics on the number of strands by length, and the number of distinct virtual helices visited by each strand.

#### Find A Specific Sequence

Source: [`examples/seq-search.py`](examples/seq-search.py)  
This script loads a Cadnano file and searches through all of the domains to find a specific DNA sequence. It will print out all domain ids that contain that sequence.

### Example Usage: Scripts

There are currently three main scripts available: 
* `converter`: Converts between different formats.
* `stapler`: Adds staples to a design following certain templates.
* `vis`: Visualize a design. Requires PyOpenGL installed (`pip install pyopengl`).

For all of these, it's assumed that you are running from the `scripts/` subdirectory.

### Running Tests: via Docker 

You can easily run the tests included in this repo. These tests should be run, and pass, before submitting an update to the repo. 
The tests are packaged in a `docker-compose` file for easy execution. Docker and Docker Compose must be installed to run tests using this method. 
See the links at the end of this section for installation links.

To run the tests, from the repo root, issue this command:
`docker-compose -f docker-compose-test.yml run test`


Here is more information on installing docker and docker-compose.
* [Installing Docker](https://docs.docker.com/engine/installation/)
* [Installing Docker Compose](https://docs.docker.com/compose/install/)

#### Converter

For more help on the converter's options:

```shell
./converter.py --help
```

To convert a Cadnano design file `my_sample.json` into a [Nanodesign Web Application](https://autode.sk/nanodesign) viewer file, using the scaffold sequence `M13mp18`:

```shell
./converter.py --infile my_sample.json --informat cadnano --inseqname M13mp18 --outfile my_sample_viewer.json --outformat viewer
```

To convert a Cadnano design file `my_sample.json` into a [CanDo](https://cando-dna-origami.org) format file, using the scaffold sequence `p8064':

```shell
./converter.py -i my_sample.json -if cadnano -isn p8064 -o my_sample.cndo -of cando
```

Note that we've abbreviated the option names here, using the short options.

Finally, to convert a Cadnano design file `my_sample.json` into a [mmCIF](http://mmcif.wwpdb.org/) file, using the sequence file `my_sample.csv` generated via Cadnano sequence export:

```shell
./converter.py -i my_sample.json -if cadnano --inseqfile my_sample.csv -o my_sample.cif -of cif
```

For using the complex options, like `--transform`, `--staples`, look at the test scripts found in `tests/converters/` for more details.

#### Stapler

The stapler currently does not have a full set of arguments or help options. It takes a Cadnano file which has a single scaffold, and already has a maximal staple set in the design file. You could generate this using the converter's staple operations, or in Cadnano. The template used is one which is optimized for honeycomb lattice structures, and the timescales are such that it should run relatively quickly, but may not be the fully optimal staple set. See the comments in `scripts/stapler.py` for how to modify these settings. 

To run the stapler on the Cadnano file `my_sample.json`:

```shell
./stapler.py my_sample.json
```

Please see the [Getting Help](#get-help) section for a forum link should you need more assistance in using the stapler!

#### Visualizer

The visualizer requires PyOpenGL installed in order to run. To install PyOpenGL:

```shell
pip install pyopengl
```

For help on the visualizer's arguments and options:

```shell
./vis.py --help
```

To visualize a Cadnano file named `my_sample.json`, using the sequence `M13mp18`:

```shell
./vis.py -i my_sample.json -isn M13mp18
```

For more details on extra options, check out the test scripts in `tests/visualizer/`.



--------

## The Package

This package has several basic components:

* Data Format - our internal data format for use by specific algorithms and to serve as an intermediate if needing to do format conversion. This format breaks a nanostructure down into core components: contiguous helices, nucleic acid strands, domains within a strand.
* Converters - we want to be able to work with any format found in the community. Currently we support loading of Cadnano files, and writing of several different formats, including [CanDo](http://cando-dna-origami.org), [SimDNA](https://github.com/mingqiu/SimDNA), [Nanodesign Viewer](https://autode.sk/nanodesign), and more.
* Algorithms - common algorithms that are used on these structures. Currently this includes melting point calculations for domains, automatic generation of staple sets for a structure, removal of staples matching certain selections, and more.

### Get help
 - API Documentation: We currently have help docstrings for most functions and methods, and will be expanding this into a full API Docs page soon.
 - [Forums](https://forum.bionano.autodesk.com/c/Nano-Design/nanodesign-python)
 - [Email us](mailto:nanodesign@autodesk.com)

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for details. 

## Contributors

We maintain an informal list of contributors to the repository. See [CONTRIBUTORS.md](CONTRIBUTORS.md) for a list of those who have made this possible.

## Notices

See [NOTICES.md](NOTICES.md) for details about incorporated code and design files.

## License

Copyright 2016 Autodesk Inc.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

