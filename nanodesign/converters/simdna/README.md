SimDNA Pairs Format

Default extension: .pairs
Encoding: ASCII text
Lines: 1st line is just an integer, indicating number of base records in the file.
2nd and remaining lines are base records, which consist of 7 columns, space separated. 

1st column: Chain (strand) id. Integer, 0-indexed
2nd column: Base id. Integer, 1-indexed. Strand-relative
3rd-5th column: X,Y,Z coordinates of C5' atom, floating point.
6th column: Chain (strand) id of paired base. Integer, 0-indexed. -1 if base is unpaired
7th column: Base id of paired base. Integer, 1-indexed. Strand-relative. -1 if base is unpaired.

IMPORTANT:  The base records MUST be sorted by chain (ascending), and then by base id within that chain (ascending). 
The chain with id 0 SHOULD be the scaffold. It's unclear how this will work with multiple scaffolds, but 0 should always be the scaffold if possible.
All bases are indexed using strand-relative base ids. So a strand with length 32 will have base ids 1-32, etc. 
