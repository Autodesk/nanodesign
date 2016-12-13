Nanodesign
==========

Create, interact and modify nucleic acid based structures. The Nanodesign Python Package is aimed at proving a toolkit for working with structural DNA/RNA nanotechnology designs. This package is meant to support all types of interaction with these designs, including loading Cadnano files and modifying them, building new structures from scratch, or just converting formats so you can easily visualize or simulate your design.

## The Package

This package have several basic components:

* Data Format - our internal data format for use by specific algorithms and to serve as an intermediate if needing to do format conversion. This format breaks a nanostructure down into core components: contiguous helices, nucleic acid strands, domains within a strand.
* Converters - we want to be able to work with any format found in the community. Currently we support loading of Cadnano files, and writing of several different formats, including [CanDo](http://cando-dna-origami.org), [SimDNA](https://github.com/mingqiu/SimDNA), [Nanodesign Viewer](https://autode.sk/nanodesign), and more.
* Algorithms - common algorithms that are used on these structures. Currently this includes melting point calculations for domains, automatic generation of staple sets for a structure, removal of staples matching certain selections, and more.

## Get help
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

