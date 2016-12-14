## Directory Structure

We use the following directory structure:

* `nanodesign`: the main python package
* `docs`: doc generation for the package
* `tests`: tests for the package
* `scripts`: example scripts for using the package, as well as scripts used to run specific tasks for the NanoDesign web application

Individual top level files of interest:

* `README.md`: You're reading it.
* `Makefile`: Some basic commands you might use, such as 'make init' for doing the requirements install, or 'make tests' for running the tests. More will be added.
* `requirements.txt`: pip requirements file. Add specific package requirements here.
* `setup.py`: distutils based package install. To run the package install, 'python setup.py install'
* `LICENSE`: Our distribution license. Currently blank, do not distribute.
* `Dockerfile`: Defines docker container with appropriate requirements.

## Directory Structure, Part 2

========================== Directory Structure ==================== 

nanodesign/

    base.py
    dna_structure.py         
    model.py
    parameters.py
    protein_structure.py
    sequence.py
    strand.py

    converters/
        converter.py

        cadnano/
            common.py
            convert_design.py ----->  NanodDnaStructure
            design.py          
            reader.py ------------->  CadnanoDesign

        cando/
            writer.py ------> .cndo 

        viewer/
            writer.py ------> .json 


========================== Converting Files ==================== 

converter.py --infile=../data/fourhelix.json --informat="cadnano"  --outfile=./results/fourhelix_viewer.json  --outformat="viewer"

converter.py --infile=../data/fourhelix.json --informat="cadnano"  --outfile=./results/fourhelix_cando.cndo   --outformat="cando"

cadnano
-------
  - read cadnano json and csv file: reader.py        
  - store data into CadnanoDesign object: design.py
  - convert data and store into NanodDnaStructure object: convert_design.py

viewer 
------
  - extract data from NanodDnaStructure object and write it to DNA Design viewer json: writer.py

cando  
-----
  - extract data from NanodDnaStructure object and write it to CanDo .cndo file: writer.py


