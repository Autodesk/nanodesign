
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


