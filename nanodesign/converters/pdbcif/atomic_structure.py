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

"""
This module defines the classes used to create an atomic representation of a DNA structure. 

The atomic representation of a DNA structure is generated from an input DnaStructure object
which contains the bases and strands making up a DNA design. An AtomicStructureStrand object
is created for each strand in the DNA structure. It is similar to the DnaStrand object but is
augmented with the rotation and translation data of each base along a strand.

A Molecule object is created for each AtomicStructureStrand object. It contains a list of
atoms for a strand.

Atomic models are generated using template structures containing three paired residues for A-T, G-C, C-G and T-A.
"""
from collections import OrderedDict
import logging
import numpy as np
import os
from math import sqrt,cos,acos,sin,asin,pi

import pdb_reader
from ...data.parameters import DnaBaseNames

# temp code to handle objects as they are being transitioned into the main package

dna_pdb_templates_dir = os.path.abspath( os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../res/'))

class Atom(object):
    """ This class stores the data for an atom.

        Attributes:
            chainID (string): The ID of the chain the atom belongs to.
            coords (NumPy 3x1 array of floats): The atom coordinates. 
            element (string): The atom element name (e.g. 'P').
            id (int): The atom ID.
            name (string): The atom name (e.g. 'OP1').
            res_name (string): The name of the residue the atom belongs to (e.g. 'DT').
            res_seq_num (int): The residue sequence number. 
    """
    def __init__( self, id, name, res_name, chainID, res_seq_num, x, y, z, element ):
        self.id = int(id)
        self.res_name = res_name
        self.res_seq_num = int(res_seq_num)
        self.chainID = chainID
        self.name = name
        self.element = element
        self.coords = np.array([float(x),float(y),float(z)],dtype=float)

    def dupe(self):
        """ Create a new Atom object with the same data. """
        return Atom(self.id, self.name, self.res_name, self.chainID, self.res_seq_num, self.coords[0], self.coords[1],
                    self.coords[2], self.element )

class MoleculeType:
    """ This class defines atomic structure types. """
    UNKNOWN = 'unknown'
    ION = 'ion'
    NUCLEIC_ACID = 'nucleic_acid'
    PROTEIN = 'protein'
    WATER = 'water'

class Molecule(object):
    """ This class stores data for a molecular structure.

        Attributes:
            atoms (List[Atom]): The atoms in this molecule.
            chains (Set[string]): The chains in this molecule.
            model_id (int): The ID of the model this molecule belongs to.
    """
    def __init__(self, model_id, type=MoleculeType.UNKNOWN):
        self.id = model_id
        self.model_id = model_id
        self.type = type 
        self.atoms = []
        self.chains = set([])
        self.residues = OrderedDict()

    def add_atom(self, atom):
        """ Add an atom to the list of atoms. 
        """
        self.atoms.append(atom)
        self.chains.add(atom.chainID)
        if atom.res_seq_num not in self.residues:
            #self.residues[atom.res_seq_num] = []
            self.residues[atom.res_seq_num] = {}
        self.residues[atom.res_seq_num][atom.name.strip()] = atom

class AtomicStructureStrand(object):
    """ This class stores data for the atomic structure of a strand. 

        Attributes:
            chainID (string): The chain ID for this strand. This is used when writing out the atomic structure.
            dna_strand (DnaStrand): The DnaStrand object representing a series of bases in the DNA structure.
            id (int): The strand ID: between 0 and number of strands-1.
            is_circular (bool): True if the strand is circular.
            is_main (List[bool]): For each base in the strand True if the base is part of a scaffold.
            rotations (List[NumPy 3x3 array of floats]): The rotation matrix for each base in the strand.
            seq (List[string]: The base sequence letter (A,C,G or T) for the strand bases.
            tour (list[int]): The list of base IDs for the strand.
            translations (List[NumPy 3x1 array of floats]): The translation vector for each base in the strand.

        This class is similar to the DnaStrand class object and uses data from that object but contains additional 
        data (e.g. rotation and translation data of each base along a strand) needed to generate an atomic structure 
        of a strand.

        The chain ID for a strand is generated using the format:

            <sc|st>.<vhelixNum>.<startPos>
        
            where sc=scaffold, st=staple
                  vhelixNUm = the virtual helix number from cadnano
                  startPos = the position in the virtual helix of the first base in the strand.
    """
    def __init__(self, strand): 
        """ Initialize a AtomicStructureStrand object. 

            Arguments:
                strand (DnaStrand): The DnaStrand object that will be used to create an atomic structure of a strand. 

        """
        self.dna_strand = strand
        self.id = strand.id
        self.tour = strand.tour
        self.is_circular = strand.is_circular
        self.is_scaffold =  strand.is_scaffold

        # Initialize rotations and translations along the strand.
        size = len(strand.tour)
        self.is_main = [False]*size
        self.seq = ['N']*size
        self.rotations = [ np.empty(shape=(0,0)) ] * size
        self.translations = [ np.empty(shape=(0)) ] * size

        # Create the chain ID.
        base_conn = strand.dna_structure.base_connectivity
        first_base = strand.tour[0]
        if strand.is_scaffold:
            self.chainID = "sc"
        else:
            self.chainID = "st"
        self.chainID = self.chainID + "." + str(first_base.h) + "." + str(first_base.p)

    def get_base_index(self, base):
        """ Get the index of the base along the strand. """
        return self.dna_strand.get_base_index(base)

class AtomicStructure(object):
    """ This class stores the atomic structure for a DNA model. 

        Attributes:
            dna_structure (DnaStructure): The DNA structure the atomic stucture will be created from. 
            molecules (List[Molecule]): The list of Molecule objects created for the atomic structure.
            strands (List[AtomicStructureStrand]): The list of atomic structure strands created for the atomic structure. 
    """ 

    # Set the names of the PDB files used as templates for generating atomic models of paired bases.
    TEMPLATE_PDB_STRUCTURE_FILE_A = 'AAA.pdb'
    TEMPLATE_PDB_STRUCTURE_FILE_G = 'GGG.pdb'
    TEMPLATE_PDB_STRUCTURE_FILE_C = 'CCC.pdb'
    TEMPLATE_PDB_STRUCTURE_FILE_T = 'TTT.pdb'

    def __init__(self, dna_structure): 
        """ Initialize a AtomicStructure object. 

            Arguments:
                dna_structure (DnaStructure): The DNA structure the atomic stucture will be created from. 
        """ 
        self.dna_structure = dna_structure 
        self.molecules = [] 
        self.strands = []
        self._logger = logging.getLogger(__name__)   
        self._init_strand_data()

        # Template structures file names.
        self.templatePdbStructureFile_A = os.path.join(dna_pdb_templates_dir, 'AAA.pdb')
        self.templatePdbStructureFile_G = os.path.join(dna_pdb_templates_dir, 'GGG.pdb')
        self.templatePdbStructureFile_C = os.path.join(dna_pdb_templates_dir, 'CCC.pdb')
        self.templatePdbStructureFile_T = os.path.join(dna_pdb_templates_dir, 'TTT.pdb')

    def _init_strand_data(self):
        """ Create the list of AtomicStructureStrand objects from DnaStrand objects. """
        for dna_strand in self.dna_structure.strands:
            atomic_strand = AtomicStructureStrand(dna_strand)
            self.strands.append(atomic_strand)

    def generate_structure(self):
        """ Generate the atomic structure for the dna model. """
        base_conn = self.dna_structure.base_connectivity
        base_nodes = self.dna_structure.helix_axis_coords
        triads = self.dna_structure.helix_axis_frames
        id_nt = self.dna_structure.id_nt
        self._logger.info("Generate atomic structure.") 
        self._logger.info("Number of bases  %d " % len(base_conn))
        self._logger.info("Number of base helix axis nodes %d " % len(base_nodes))
        self._logger.info("Number of triads: %d" % triads.shape[2])

        # Convert base node coords to angstroms.
        for i in xrange(0,len(id_nt)):
            base_nodes[i,0] = 10.0*base_nodes[i,0] 
            base_nodes[i,1] = 10.0*base_nodes[i,1] 
            base_nodes[i,2] = 10.0*base_nodes[i,2] 

        # Reverse the directions of coordinate frames for the major groove (0) and helix axis (2). 
        if True:
            for i in xrange(0,len(id_nt)):
                triads[0,0,i] = -triads[0,0,i] 
                triads[1,0,i] = -triads[1,0,i] 
                triads[2,0,i] = -triads[2,0,i]
                triads[0,2,i] = -triads[0,2,i]
                triads[1,2,i] = -triads[1,2,i]
                triads[2,2,i] = -triads[2,2,i]

        # Rotation to convert frame [e1 e2 e3] to [-e2 e3 -e1]. 
        rot_mat = np.array([ [0, 0, -1], [ -1, 0, 0], [0, 1, 0]], dtype=float)

        # Set the rotation matrix, translation vector and sequence for strand paired bases.
        self._logger.debug("=================== set strand data ==================");
        for i in xrange(0,len(id_nt)):
            self._logger.debug("---------- i %d ----------\n" % i)
            base_id1 = id_nt[i,0] 
            base1 = base_conn[base_id1]
            strand1 = self.strands[base1.strand]
            base1_sindex = strand1.get_base_index(base1)
            strand1.is_main[base1_sindex] = True;
            strand1.seq[base1_sindex] = base1.seq
            strand1.rotations[base1_sindex] = np.dot(triads[:,:,i],rot_mat)
            strand1.translations[base1_sindex] = base_nodes[i]
            #self._logger.debug(">>> nodes  %s" % str(base_nodes[i])) 
            #self._logger.debug(">>> triad 1: %6.2g %6.2g %6.2g" % (triads[0,0,i], triads[1,0,i], triads[2,0,i])) 
            #self._logger.debug("          2: %6.2g %6.2g %6.2g" % (triads[0,1,i], triads[1,1,i], triads[2,1,i])) 
            #self._logger.debug("          3: %6.2g %6.2g %6.2g" % (triads[0,2,i], triads[1,2,i], triads[2,2,i])) 
            #sys.exit(0)

            # Paired base.
            base_id2 = id_nt[i,1]             
            base2 = base_conn[base_id2]
            strand2 = self.strands[base2.strand]
            base2_sindex = strand2.get_base_index(base2)
            strand2.is_main[base2_sindex] = False;
            strand2.seq[base2_sindex] = base2.seq
            strand2.rotations[base2_sindex] = np.dot(triads[:,:,i],rot_mat)
            strand2.translations[base2_sindex] = base_nodes[i]
        #__for i in xrange(0,len(id_nt1))

        self._logger.debug("=================== create bulges ==================");
        for strand in self.strands:
            self._logger.debug("=================== strand:%d =================== "  % strand.id)
            base = base_conn[strand.tour[0]-1]
            rotations = strand.rotations
            translations = strand.translations
            is_main = strand.is_main
            is_circular = strand.is_circular
            seq = strand.seq

            self._generate_bulge_dof(strand, rotations, translations, is_main, is_circular)

            strand.rotations = rotations 
            strand.translations = translations
            strand.is_main = is_main

            for i in xrange(0,len(seq)):
                if (seq[i] == 'N'):
                    base_index = strand.tour[i]-1
                    strand.seq[i] = base_conn[base_index].seq
            #_for i in xrange(0,len(seq))
        #__for strand in self.strands

        # Generate atomic structures from the dna strands.
        self._logger.debug("=================== generate atomic structures ==================");
        self.molecules = self._generate_atoms(self.strands)
        self._logger.debug("Generated %d atomic structures. " % len(self.molecules));
        return self.molecules 

    def _generate_atoms(self, strands):
        """ Generate atomic structures from the dna strands. 

            Arguments:
                strands (List[AtomicStructureStrand]): The list of atomic stucture strands to create atoms for.
        """
        # Read template structures, seperating atoms into forward (5'->3') and reverse chains.
        A_for, T_rev = self._read_template(AtomicStructure.TEMPLATE_PDB_STRUCTURE_FILE_A)
        G_for, C_rev = self._read_template(AtomicStructure.TEMPLATE_PDB_STRUCTURE_FILE_G)
        C_for, G_rev = self._read_template(AtomicStructure.TEMPLATE_PDB_STRUCTURE_FILE_C)
        T_for, A_rev = self._read_template(AtomicStructure.TEMPLATE_PDB_STRUCTURE_FILE_T)

        # Create a dict mapping base name to forward and reverse structures. 
        forward_struct = { DnaBaseNames.A : A_for, DnaBaseNames.C : C_for, DnaBaseNames.G : G_for, DnaBaseNames.T : T_for}
        reverse_struct = { DnaBaseNames.A : A_rev, DnaBaseNames.C : C_rev, DnaBaseNames.G : G_rev, DnaBaseNames.T : T_rev}

        molecular_structures = []
        num_strands = len(strands)
        first_atomID = 1     # serial number for the first atom in a tour

        # Create a Molecule object containing the atoms for each strand.
        for strand in strands:
            molecule, first_atomID = self._create_atoms_from_strand(strand, forward_struct, reverse_struct, first_atomID)
            molecular_structures.append(molecule)
        #__for strand in strands

        return molecular_structures 
    #__def _pdb_generate(self, strands)

    def _create_atoms_from_strand(self, strand, forward_struct, reverse_struct, first_atomID_in):
        """ Create the atoms for the list of bases in a strand.

            Arguments:
                strand (AtomicStructureStrand): The input strand to create and atomistic structure from. 
                forward_struct: 
                reverse_struct: 
                first_atomID_in:

            Returns: 
                molecule (Molecule):  A Molecule object containing the strand atoms.

            A Molecule object is created to store the atoms for the strand structure. 
        """
        self._logger.debug("=================== strand ID %d =================== " % strand.id ) 
        self._logger.debug("Number of bases %d " % len(strand.tour) ) 
        self._logger.debug("Scaffold %d " % strand.is_scaffold ) 
        base_conn = strand.dna_strand.dna_structure.base_connectivity
        base = base_conn[strand.tour[0]-1]
        self._logger.debug( "Start helix %d  position %d" % (base.h, base.p) ) 
        first_atomID_out = first_atomID_in
        max_atom_per_base = 40
        num_bases = len(strand.tour)

        # Create a Molecule object to store the atoms for the strand structure. 
        molecule = Molecule(strand.id)

        # Create a rotation matrix for a 180 degree rotation about the y-axis.
        Ry180 = self._Ry(180.0)

        # Build up a structure from the strand bases.
        atom_index = 1
        for i in xrange(0,num_bases):
            strand_R = strand.rotations[i]
            if (strand_R.size == 0):
                continue
            base_name = strand.seq[i].upper()
            self._logger.debug( "base(%d)='%s'" % (i,base_name)) 

            if base_name not in forward_struct:
                self._logger.warn( "base(%d)='%s' not found." % (i,base_name)) 
                continue 

            # Get the rotation and translation for the base.
            R = np.dot(strand_R,Ry180)
            D = strand.translations[i]

            # Get the structure for the base. If is_main[i] is True then the base is for a scaffold.
            if (strand.is_main[i]):       
                current_struct = forward_struct[base_name]    
            else: 
                current_struct = reverse_struct[base_name]   

            # Duplicate and transform the atoms in the reference structure.
            num_atoms = len(current_struct)
            for j in xrange(0,num_atoms):
                current_atom = current_struct[j].dupe()
                current_atom.id = first_atomID_out
                current_atom.chainID = strand.chainID
                first_atomID_out = first_atomID_out+1
                current_atom.res_seq_num = i+1;
                coords = current_atom.coords
                # Transform the atom coordinates. 
                xform_coord = current_atom.coords
                xform_coord = np.dot(R,xform_coord)
                xform_coord = np.add(xform_coord, D)
                self._logger.debug(" atom id %d  name %s  coord (%g %g %g)" % (current_atom.id, current_atom.name,
                    xform_coord[0], xform_coord[1], xform_coord[2]))
                current_atom.coords = xform_coord
                molecule.add_atom(current_atom)
                atom_index = atom_index+1
            #__for j in xrange(0,num_atoms):
        #_for i in xrange(0,num_bases)

        return molecule, first_atomID_out 

    def _read_template(self, infile_name):
        """ Read a template structure from a PDB file.

            Arguments:
                infile_name (String): The name of the template structure file to read.

            Returns: 
                forward_struct (List[Atom]): The list of atoms for the forward dna strand.
                reverse_struct (List[Atom]): The list of atoms for the reverse dna strand.

            In the reference structure:
                forward strand chainID = S
                backward strand chainID = A
        """ 
        forward_struct = []
        reverse_struct = []

        # Read in the template structure.
        file_name = os.path.join(dna_pdb_templates_dir, infile_name) 
        reader = pdb_reader.PdbReader()
        reader.read(file_name)
        molecules = reader.molecules
        molecule = molecules[0]

        # Create rotation matrix for a -90 deg rotation about x-axis
        # followed by a -90 deg rotation about y-axis.
        R = np.dot(self._Ry(-90), self._Rx(-90))

        # Group atoms into forward/backward chains.
        for atom in molecule.atoms: 
            # Transform and store atom for the forward reference structure.
            if (atom.chainID == 'S') and (atom.res_seq_num == 2):
                atom.id = -1;
                atom.res_seq_num = -1;
                xform_coord = atom.coords
                xform_coord = np.dot(R, xform_coord)
                atom.coords = xform_coord
                atom.chainID = 'A'
                forward_struct.append(atom)

            # Transform and store atom for the backward reference structure.
            elif (atom.chainID == 'A') and (atom.res_seq_num == 2):
                atom.id = -1;
                atom.res_seq_num = -1;
                xform_coord = atom.coords
                xform_coord = np.dot(R,xform_coord)
                atom.coords = xform_coord
                reverse_struct.append(atom)
        #__for atom in molecule.atoms: 

        return forward_struct, reverse_struct 

    def _generate_bulge_dof(self, strand, rotations, translations, is_main, is_circular):
        self._logger.debug(">>> _generate_bulge_dof ")
        num_nt = len(rotations)

        # Find all regions in the strand that are not paired.
        single_strands = self._find_single_strands(rotations, is_circular)

        # Remove single-stranded regions with free ends.
        if (not is_circular):
            for i in xrange(len(single_strands)-1,-1,-1):
                strand_pos = single_strands[i]
                if ( (strand_pos[0] == 0) or (strand_pos[-1] == num_nt-1) ):
                    single_strands[i] = [];
        #_if (not is_circular)

        self._logger.debug("Number of single_strands) %d " % len(single_strands))

        # Find positions and orientations of the single-stranded regions.
        for i in xrange(0,len(single_strands)):
            strand_pos = single_strands[i]
            if len(strand_pos) == 0:
                continue 
            self._logger.debug( ">>> process single stranded region %d " % i) 
            self._logger.debug( "strand_pos %d: %s " % ( len(strand_pos), str(strand_pos)))
            #sys.exit(0)

            if self._logger.getEffectiveLevel() == logging.DEBUG:
                base_conn = strand.dna_strand.dna_structure.base_connectivity
                base_str = ""
                for si in strand_pos:
                    base_id = strand.tour[si]-1
                    base = base_conn[base_id]
                    base_str += "(%d,%d) " % (base.h,base.p)
                self._logger.debug( "bases: %s" % base_str)
            #__if self._logger.getEffectiveLevel() == logging.DEBUG

            # Get the rotation/translation at the start of the single-strand region.
            i_1 = strand_pos[0] - 1
            if (i_1 < 0):
                i_1 = num_nt-1
            R_1 = rotations[i_1]
            d_1 = translations[i_1]

            # Get the rotation/translation at the enb of the single-strand region.
            if (not is_main[i_1]):
                vrot = self._vrrotvec2mat([0, 0, 1], pi)
                R_1 = np.dot(R_1,vrot)
            i_2 = strand_pos[-1] + 1
            if (i_2 > num_nt-1):
                i_2 = 0
            R_2 = rotations[i_2]
            d_2 = translations[i_2]
            self._logger.debug( "i_1 %d   i_2 %d " % (i_1, i_2)) 

            if (not is_main[i_2]):
                vrot = self._vrrotvec2mat([0, 0, 1], pi)
                R_2 = np.dot(R_2,vrot)
    
            # Interpolate between the start/end rotation/translation.
            R_fit, d_fit = self._fit_R_d(R_1, d_1, R_2, d_2, len(strand_pos))
    
            # Update the rotations and translations.
            for j in xrange(0,len(strand_pos)):
                k = strand_pos[j]
                is_main[k] = True
                rotations[k] = R_fit[:,:,j]
                translations[k] = d_fit[:,j]
            #__for j in xrange(0,len(strand_pos))

        #__for i in xrange(0,len(single_strands))

    #__def _generate_bulge_dof

    def _find_single_strands(self, rotations, is_circular):
        """ Find all regions in the strand that are not paired.

            Arguments:
                rotations (List[NumPy 3x3 array of floats]): The rotation matrix for each base in the strand.
                is_circular (bool): If True then the strand is circular.

            Returns: 
                single_strands (List(List[2])): A list of [start,end] single-strand locations within the strand.

            Iterate through the strand 'rotations' array  to find bases that are not paired.
        """
        num_nt = len(rotations)
        is_visited = [False]*num_nt

        for i in xrange(0,num_nt):
            rmat = rotations[i]
            if (rmat.size != 0):
                is_visited[i] = True;
        #_for i

        # Search for the start of a single-stranded region.
        loop = 0
        single_strands = []

        while (True):
            try:
                i_0 = is_visited.index(False)
            except:
                break

            single_strands.append([])

            # Scan in the 5' direction. 
            i = i_0
            while (not is_visited[i]):
                is_visited[i] = True;
                single_strands[loop].append(i);

                if (i > 0):
                    i = i - 1
                elif (is_circular):
                    i = num_nt-1
                else:
                    break
            #_while(not is_visited[i])

            single_strands[-1] = []
            is_visited[i_0] = False

            # Scan in the 3' direction.
            i = i_0
            while(not is_visited[i]):
                is_visited[i] = True
                single_strands[loop].append(i)

                if (i < num_nt-1):
                    i = i + 1
                elif (is_circular):
                    i = 0
                else:
                    break;
            #_while(not is_visited[i])

            loop = loop + 1
        #_while(True)

        return single_strands
    #__def _find_single_strands

    def _vrrotvec2mat(self, axis, theta):
        """ Create a rotation matrix to rotate theta degrees about the axis defined by vec. """
        s = np.sin(theta)
        c = np.cos(theta)
        t = 1 - c
        #print("[vrrotvec2mat] theta=%s" % str(theta))

        # normalize the vector
        x,y,z = axis / np.linalg.norm(axis)
        #print("[vrrotvec2mat] x=%s" % str(x))
        #print("[vrrotvec2mat] y=%s" % str(y))
        #print("[vrrotvec2mat] z=%s" % str(z))

        return np.array([ [t*x*x + c,   t*x*y - s*z,  t*x*z + s*y],
                          [t*x*y + s*z, t*y*y + c,    t*y*z - s*x],
                          [t*x*z - s*y, t*y*z + s*x,  t*z*z + c  ]])

    def _fit_R_d(self, R_1, d_1, R_2, d_2, num_bases):

        """ This function interpolates the rotations and translations between the start and end of a region in a strand. 

            Arguments:
                R_1 (NumPy 3x3 array of floats): The rotation matix at the start of the strand region. 
                d_1 (NumPy 3x1 array of float): The translation vector at the start of the strand region.
                R_2 (NumPy 3x3 array of floats): The rotation matix at the end of the strand region. 
                d_2 (NumPy 3x1 array of float): The translation vector at the end of the strand region.
                num_bases (int): The number of bases in the strand region.

            Returns: 
                R_fit (List[NumPy 3x3 array of floats)]: The list of rotation matices from the start to the end 
                    of the strand region. 
                d_fit (List[NumPy 3x1 array of floats)]: The list of translation matices from the start to the end 
                    of the strand region. 

            Given the rotations and translations at the start of a strand region (R_1,d_1) and the rotations and translations 
            at the end of a strand region (R_2,d_2) we calculate the intermediate (R,d) for each base in between. To interpolate
            rotations we first calculate the matrix R that rotates R_1 into R_2

                R * R_1 = R_2 
                R = R_2 * transpose(R_1)

            The angle of the (axis,angle) representation of R is then used to generate successive rotation matices. 

        """ 
        R_fit = np.zeros((3,3,num_bases), dtype=float)
        d_fit = np.zeros((3,num_bases), dtype=float)

        # Calculate R from right-inverse of R * R_1 = R_2. 
        R = np.dot(R_2,R_1.T)

        # Extract angle/axis of rotation.
        axis1,theta1 = self._vrrotmat2vec(R_1);
        axis2,theta2 = self._vrrotmat2vec(R_2);
        axis,theta = self._vrrotmat2vec(R);
        self._logger.debug( "[_fit_R_d] theta  %g  axis  %g %g %g" % (theta, axis[0], axis[1], axis[2]))
        self._logger.debug( "[_fit_R_d] theta1 %g  axis1 %g %g %g" % (theta1, axis1[0], axis1[1], axis1[2]))
        self._logger.debug( "[_fit_R_d] theta2 %g  axis2 %g %g %g" % (theta2, axis2[0], axis2[1], axis2[2]))

        # Interpolate rotations/translations. 
        for i in xrange(0,num_bases):
            j = i+1
            angle = (theta*j)/(num_bases+1)
            self._logger.debug( "[_fit_R_d] angle %g " % angle) 
            R = self._vrrotvec2mat(axis, angle)
            R_fit[:,:,i] = np.dot(R,R_1)
            d_fit[:,i] = (d_1*(num_bases+1-j) + d_2*j) / (num_bases+1)
        #__for i in xrange(0,num_bases)

        return R_fit, d_fit

    def _vrrotmat2vec(self,R):
        """ Extract the equivalent rotation about an axis from a rotation matrix. """
        #print("-------------- _vrrotmat2vec -------------- ")
        r00 = R[0,0]
        r01 = R[0,1]
        r02 = R[0,2]
        r10 = R[1,0]
        r11 = R[1,1]
        r12 = R[1,2]
        r20 = R[2,0]
        r21 = R[2,1]
        r22 = R[2,2]
        #print("R %s " % str(R))
        #print("r21-r12 %s " % str(r21-r12))
        #print("r02-r20 %s " % str(r02-r20))
        #print("r10-r01 %s " % str(r10-r01))

        angle = acos(( r00 + r11 + r22 - 1)/2.0)
        #print("angle %g " % angle) 

        # Check for singularities at angle=0.
        if abs(angle) < 0.0001: 
            angle = 0.0
            x = 1.0
            y = 0.0
            z = 0.0

        # Check for singularity at angle=180. Calculate axis using largest diagonal term.
        elif abs(angle - pi) < 0.0001: 
            epsilon = 0.01
            xx = (r00+1) / 2.0
	    yy = (r11+1) / 2.0
	    zz = (r22+1) / 2.0
	    xy = (r01+r10) / 4.0
	    xz = (r02+r20) / 4.0
	    yz = (r12+r21) / 4.0

	    if (xx > yy) and (xx > zz): 
                if xx < epsilon:
                    x = 0.0
		    y = 0.7071
		    z = 0.7071
                else:
		    x = sqrt(xx)
		    y = xy / x
		    z = xz / x
	    elif yy > zz:
                if yy < epsilon: 
		    x = 0.7071
	            y = 0.0
	            z = 0.7071
		else:
		    y =sqrt(yy)
		    x = xy / y
		    z = yz / y
            else:
                if zz < epsilon: 
                    x = 0.7071
		    y = 0.7071
		    z = 0.0
                else:
		    z = sqrt(zz)
		    x = xz / z
		    y = yz / z

        # Calculate axis for non-singular angle.
        else:
            x = (r21 - r12) / sqrt( pow(r21-r12,2) + pow(r02-r20,2) + pow(r10-r01,2) )
            y = (r02 - r20) / sqrt( pow(r21-r12,2) + pow(r02-r20,2) + pow(r10-r01,2) )
            z = (r10 - r01) / sqrt( pow(r21-r12,2) + pow(r02-r20,2) + pow(r10-r01,2) )

        # Check the result.
        check_result = False
        check_result = True
        if check_result:
            c = cos(angle)
            s = sin(angle)
            t = 1.0 - c
            M = np.array([ [0, 0,-1], [ -1, 0, 0], [0, 1, 0]], dtype=float)
            M[0,0]= t*x*x + c;        M[0,1] = t*x*y - z*s;       M[0,2] = t*x*z + y*s;
            M[1,0]= t*x*y + z*s;      M[1,1] = t*y*y + c;         M[1,2] = t*y*z - x*s;
            M[2,0]= t*x*z - y*s;      M[2,1] = t*y*z + x*s;       M[2,2] = t*z*z + c;

            for i in xrange(0,3):
               for j in xrange(0,3):
                   d = R[i,j] - M[i,j]
                   if abs(d) > 0.0001:
                       print("[_vrrotmat2vec] **** WARNING: R-M is not zero: %g " % (R[i,j]-M[i,j]))
            #__for i in xrange(0,3)
        #__if check_result

        return np.array([x,y,z],dtype=float),angle

    def _Rx(self, deg):
        """ Generate a rotation matrix for an angle about the X-axis. """
        rad = pi*(deg /180.0)
        c = cos(rad)
        s = sin(rad)
        R = np.array( [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]] )
        return R

    def _Ry(self, deg):
        """ Generate a rotation matrix for an angle about the Y-axis. """
        rad = pi*(deg /180.0)
        c = cos(rad)
        s = sin(rad)
        R = np.array( [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]] )
        return R

    def get_extent(self):
        """ Get the extent of the atom coordinates. 
        """
        xmin = 0.0
        xmax = 0.0 
        ymin = 0.0
        ymax = 0.0 
        zmin = 0.0
        zmax = 0.0 
        num_atoms = 0

        if len(self.molecules) == 0:
           self.generate_structure()

        for molecule in self.molecules:
            for atom in molecule.atoms:
                x = atom.coords[0]
                y = atom.coords[1]
                z = atom.coords[2]
                if num_atoms == 0:
                    xmin = xmax = x
                    ymin = ymax = y
                    zmin = zmax = z
                else:
                    if x < xmin: xmin = x
                    elif x > xmax: xmax = x
                    if y < ymin: ymin = y
                    elif y > ymax: ymax = y
                    if z < zmin: zmin = z
                    elif z > zmax: zmax = z
                num_atoms += 1

            #__for atom in molecule.atoms
        #__for molecule in self.molecules
        return xmin,xmax,ymin,ymax,zmin,zmax 

    #========================================================================================================
    #======================================= new ssDNA generation code ======================================
    #========================================================================================================

    def generate_structure_ss(self):
        """ Generate the atomic structure for the dna model. """
        base_conn = self.dna_structure.base_connectivity
        self._logger.info("Generate atomic structure for ssDNA.") 
        self._logger.info("Number of bases  %d " % len(base_conn))

        # Scale to convert base node coords in nm to angstroms.
        nm_to_ang = 10.0

        # Rotation to convert frame [e1 e2 e3] to [-e2 e3 -e1].
        rot_mat = np.array([ [0, 0, -1], [ -1, 0, 0], [0, 1, 0]], dtype=float)

        # Create strands helix maps. 
        self.dna_structure.set_strand_helix_references()

        # Set the rotation matrix, translation vector and sequence for strand paired bases.
        helix_map = self.dna_structure.structure_helices_map
        for strand in self.strands:
            self._logger.debug("=================== strand:%d =================== "  % strand.id)
            base_coords = strand.dna_strand.get_base_coords()

            for i,base in enumerate(strand.tour):
                #frame = np.zeros((3,3),dtype=float)
                #frame = np.copy(base.ref_frame)
                frame = base.ref_frame
                self._logger.debug("base id %d  vh %d  pos %d " % (base.id, base.h, base.p))
                frame[:,0] = -frame[:,0]
                frame[:,2] = -frame[:,2]
                self._logger.debug("     frame[0] %g %g %g" % (frame[0,0], frame[1,0], frame[2,0])) 
                self._logger.debug("     frame[1] %g %g %g" % (frame[0,1], frame[1,1], frame[2,1])) 
                self._logger.debug("     frame[2] %g %g %g" % (frame[0,2], frame[1,2], frame[2,2])) 
                strand.is_main[i] = strand.is_scaffold
                if base.seq == "N":
                    strand.seq[i] = "A"
                else:
                    strand.seq[i] = base.seq
                strand.rotations[i] = np.dot(frame,rot_mat)
                strand.translations[i] = nm_to_ang*base_coords[i]
            #__for base_id in strand.tour
        #__for strand in self.strands

        self._logger.debug("=================== create bulges ==================");
        for strand in self.strands:
            self._logger.debug("---------- strand %d ---------- "  % strand.id)
            base = strand.tour[0]
            rotations = strand.rotations
            translations = strand.translations
            is_main = strand.is_main
            is_circular = strand.is_circular
            seq = strand.seq

            self._generate_bulge_dof(strand, rotations, translations, is_main, is_circular)

            strand.rotations = rotations 
            strand.translations = translations
            strand.is_main = is_main

            for i in xrange(0,len(seq)):
                if (seq[i] == 'N'):
                    base_index = strand.tour[i]-1
                    strand.seq[i] = base_conn[base_index].seq
            #_for i in xrange(0,len(seq))
        #__for strand in self.strands

        # Generate atomic structures from the dna strands.
        self.molecules = self._generate_atoms_ss(self.strands)
        self._logger.debug("Generated %d atomic structures. " % len(self.molecules));
        return self.molecules 

    def _generate_atoms_ss(self, strands):
        """ Generate atomic structures from the dna strands. 

            Arguments:
                strands (List[AtomicStructureStrand]): The list of atomic stucture strands to create atoms for.
        """
        self._logger.debug("=================== _generate_atoms_ss ==================");
        # Read template structures, seperating atoms into forward (5'->3') and reverse chains.
        A_for, T_rev = self._read_template(AtomicStructure.TEMPLATE_PDB_STRUCTURE_FILE_A)
        G_for, C_rev = self._read_template(AtomicStructure.TEMPLATE_PDB_STRUCTURE_FILE_G)
        C_for, G_rev = self._read_template(AtomicStructure.TEMPLATE_PDB_STRUCTURE_FILE_C)
        T_for, A_rev = self._read_template(AtomicStructure.TEMPLATE_PDB_STRUCTURE_FILE_T)

        # Create a dict mapping base name to forward and reverse structures. 
        forward_struct = { DnaBaseNames.A : A_for, DnaBaseNames.C : C_for, DnaBaseNames.G : G_for, DnaBaseNames.T : T_for}
        reverse_struct = { DnaBaseNames.A : A_rev, DnaBaseNames.C : C_rev, DnaBaseNames.G : G_rev, DnaBaseNames.T : T_rev}

        molecular_structures = []
        num_strands = len(strands)
        first_atomID = 1     # serial number for the first atom in a tour

        # Create a Molecule object containing the atoms for each strand.
        for strand in strands:
            molecule, first_atomID = self._create_atoms_from_strand_ss(strand, forward_struct, reverse_struct, first_atomID)
            molecular_structures.append(molecule)
        #__for strand in strands

        return molecular_structures 
    #__def _pdb_generate_ss(self, strands)

    def _create_atoms_from_strand_ss(self, strand, forward_struct, reverse_struct, first_atomID_in):
        """ Create the atoms for the list of bases in a strand.

            Arguments:
                strand (AtomicStructureStrand): The input strand to create an atomistic structure from. 
                forward_struct (Dict{String:List[Atom]}: A dictionary mapping a base name to a forward structure.
                reverse_struct (Dict{String:List[Atom]}: A dictionary mapping a base name to a reverse structure.
                first_atomID_in (int): The current atom ID.

            Returns: 
                molecule (Molecule):  A Molecule object containing the strand atoms.
                first_atomID_out (int): The current atom ID of the last atom created for the input the strand. 

            A Molecule object is created to store the atoms for the strand structure. 
        """
        self._logger.debug("=================== _create_atoms_from_strand_ss =================== " ) 
        self._logger.debug("Strand %d " % strand.id ) 
        self._logger.debug("Number of bases %d " % len(strand.tour) ) 
        self._logger.debug("Scaffold %d " % strand.is_scaffold ) 
        base_conn = strand.dna_strand.dna_structure.base_connectivity
        base = strand.tour[0]
        self._logger.debug( "Start helix %d  position %d" % (base.h, base.p) ) 
        first_atomID_out = first_atomID_in
        max_atom_per_base = 40
        num_bases = len(strand.tour)
        self._logger.debug("first_atomID_in %d " % first_atomID_in) 

        # Create a Molecule object to store the atoms for the strand structure. 
        molecule = Molecule(strand.id)

        # Create a rotation matrix for a 180 degree rotation about the y-axis.
        Ry180 = self._Ry(180.0)

        # Build up a structure from the strand bases.
        atom_index = 1
        for i in xrange(0,num_bases):
            strand_R = strand.rotations[i]
            if (strand_R.size == 0):
                continue
            base_name = strand.seq[i].upper()
            self._logger.debug("---------- base(%d)='%s' ---------- " % (i,base_name)) 

            if base_name not in forward_struct:
                self._logger.warn( "base(%d)='%s' not found." % (i,base_name)) 
                continue 

            # Get the rotation and translation for the base.
            R = np.dot(strand_R,Ry180)
            D = strand.translations[i]

            # Get the structure for the base. If is_main[i] is True then the base is for a scaffold.
            if (strand.is_main[i]):       
                current_struct = forward_struct[base_name]    
            else: 
                current_struct = reverse_struct[base_name]   

            # Duplicate and transform the atoms in the reference structure.
            num_atoms = len(current_struct)
            for j in xrange(0,num_atoms):
                current_atom = current_struct[j].dupe()
                current_atom.id = first_atomID_out
                current_atom.chainID = strand.chainID
                first_atomID_out = first_atomID_out+1
                current_atom.res_seq_num = i+1;
                coords = current_atom.coords
                # Transform the atom coordinates. 
                xform_coord = current_atom.coords
                xform_coord = np.dot(R,xform_coord)
                xform_coord = np.add(xform_coord, D)
                self._logger.debug(" atom id %d  name %s  coord (%g %g %g)" % (current_atom.id, current_atom.name,
                    xform_coord[0], xform_coord[1], xform_coord[2]))
                current_atom.coords = xform_coord
                molecule.add_atom(current_atom)
                atom_index = atom_index+1
            #__for j in xrange(0,num_atoms):
        #_for i in xrange(0,num_bases)

        return molecule, first_atomID_out 

#__class AtomicStructure


