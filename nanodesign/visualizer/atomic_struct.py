#!/usr/bin/env python
""" This module is used to visualize the atomic structure of a single strand (ssDNA) of a DNA model. 

    An atomistic model of a DNA nanodesign consists of one or more scaffold strands and many staple
    strands. Each strand of DNA is represnted as a polymer chain with a unique chain ID. The VisAtomicStructure
    class in this module creates and manages the visualization of a single strand. A DNA strand is visualized 
    using several representations: 

        bonds - The atomic bonds are displayed as lines. The sugar and base planes are displayed as filled polygons.
        backbone P atoms - The DNA P atoms are displayed as spheres connected by lines.
        P atom distances - Large deviations of P-P bond lengths are highlighted.  
"""
import logging 
import numpy as np
from .geometry import VisGeometryPath,VisGeometryLines,VisGeometryPolygon
from math import sqrt

class VisAtomicStructureRepType:
    """ This class defines the atomic structure visualization representation types. """
    NONE       = 'none'
    BACKBONE   = 'backbone'
    BONDS      = 'bonds'
    CHECK      = 'check'

class VisAtomicStructure(object):
    """ This class is used to visualize the atomic structure of a DNA strand. 

        Attributes:
            color (List[float]): The color assigned to the structure. This is a list of four floats (RGBA) 
                obtained from the DNA strand color.
            graphics (VisGraphics): The VisGraphics object.
            id (int): The structure id, usually just the molecule count. Not really used.
            molecule (Molecule): The Molecule object storing the strand atomic structure.
            name (String): The strand chain name.
            representations (Dict[List[VisGeometry]): The dictionary storing the list of geometry for a representation.
            scale (Float): The scale to convert angstroms to nanometers.
            strand (VisStrand): The VisStrand object associated with this atomic structure. The color and name of the
                structure is obtained from this strand.
            strand_name (String): The strand name created for this structure.

        The visualization geometry is createed from the atomic coordinates of the atomistic model of a single DNA strand 
        taken from the given Molecule object. It is assumed that the Molecule object contains a single continuous strand
        when creating geometry.

        The chain ID of the atomic structure DNA strand is assumed to have the following format:

            <sc|st>.<vhelixNum>.<startPos>

            where sc=scaffold, st=staple
                  vhelixNUm = the virtual helix number from cadnano
                  startPos = the position in the virtual helix of the first base in the strand.

        To chain ID is modified to match the associated visualization strand name by:
            - replacing 'sc' by 'Scaffold' 
            - replacing 'st' by 'staple' 
            - replacing '.' by '_'

    """
    average_p_bond_length = 0.65   # The average distance between P-P bonds.
    p_bond_length_tol = 0.1        # The tolerance for P-P bond deviation.

    def __init__(self, id, model, molecule, graphics):
        """ Initialize a VisAtomicStructure object. 

            Arguments:
                id (int): The atomic structure ID, from 0 to the number of atomic structures.
                model (VisModel): The visualization model object used to interface with the DNA design structure
                    and manage the visualization of all representions.
                molecule (Molecule): The atomic structure molecule object for a ssDNA. 
                graphics (VisGraphics): The visualization graphics object that manages the display of geometry for
                    various representations.
        """
        self.id = id
        self.graphics = graphics
        self.molecule = molecule
        self.scale = 0.1
        self.representations = {}
        # Set the structure name. This will just be the chain ID. There should be just one chain.
        chains = list(molecule.chains)
        self.name = chains[0]
        self._logger = self._setup_logging()
        # Set the equivalent strand name for this chain.
        tokens = self.name.split(".")
        if tokens[0] == "sc":
            tokens[0] = "Scaffold"
        else:
            tokens[0] = "staple"
        self.strand_name = "_".join(tokens)
        if self.strand_name in model.strands:
            self.strand = model.strands[self.strand_name]
            self.color = self.strand.color[:]
        else:
            self.color = [1.0,1.0,1.0,1.0]
        # Set the methods to create geometry for the different representations.
        self.create_rep_methods = { 
            VisAtomicStructureRepType.BACKBONE : self.create_backbone_rep,
            VisAtomicStructureRepType.BONDS    : self.create_bonds_rep,
            VisAtomicStructureRepType.CHECK    : self.create_check_rep 
        }

    def _setup_logging(self):
        """ Set up logging."""
        logger = logging.getLogger(__name__ + ':' + self.name) 
        logger.setLevel(logging.INFO)
        # Create console handler and set format.
        console_handler = logging.StreamHandler()
        formatter = logging.Formatter('[%(name)s] %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        return logger

    @staticmethod
    def compare(a,b):
        """ The compare function used to sort strands by helix number and then position. """
        return cmp(a.strand.start_helix,b.strand.start_helix) or cmp(a.strand.start_pos,b.strand.start_pos)

    def show(self, rep, show, display=True):
        """ Show the strand atomic structure with the given representation. 

            Arguments:
                rep (String): The representation name from VisAtomicStructureRepType.
                show (bool): If true then show the geometry for the representation, else hide it.
                dispay (bool): If true then redisplay all the graphics geometry.

            If the geometry for the representation has not been created then create and store it. 
        """ 
        self._logger.debug("Show atomic structure \'%s\'  rep \'%s\' " % (self.name, rep))
        if rep not in self.representations:
            self.create_rep_methods[rep]()
        for geom in self.representations[rep]:
            geom.visible = show
        if display:
            self.graphics.display()

    def print_info(self):
        """ Print atomic stucture information. """ 
        self._logger.info("Number of residues %d" % (len(self.molecule.residues)))

    def create_backbone_rep(self):
        """ Create the geometry for the atomic structure backbone representation. 

            The strand backbone is visualized by displaying its P atoms as solid spheres connected by lines.
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("Create atomic structure backbone rep.")
        s = self.scale
        # Extract the P atom coordinates.
        points = []
        for atom in self.molecule.atoms:
            if atom.element.strip() == 'P':
                point = [s*atom.coords[0], s*atom.coords[1], s*atom.coords[2]]
                points.append(point)
        #__for atom in self.molecule.atoms
        # Create the backbone geometry.
        show_vertices = True
        name = "Backbone:%s" %  self.strand_name
        geom = VisGeometryPath(name, points, show_vertices)
        geom.line_width = 3.0
        geom.color = self.color[:]
        geom.select_vertex = True
        geom.selected_callback = self.select_backbone
        geom.entity_indexes = range(0,len(points))
        self.representations[VisAtomicStructureRepType.BACKBONE] = [geom]
        self.graphics.add_render_geometry(geom)

    def select_backbone(self, geom, index):
        """ Process backbone selection. 

            Arguments:
                geom (VisGeometry): The geometry selected.
                index (int): The index into the geometry selected.
        """
        num_res = len(self.molecule.residues)
        res_list = list(self.molecule.residues.values())
        if (index >= 0) and (index < num_res):
            residue_atoms = res_list[index]
            atom = residue_atoms["P"]
            self._logger.info("Selected backbone %s. Residue sequence number %d  Base name %s" % (self.name, 
                atom.res_seq_num, atom.res_name ))
            self.print_info()

    def create_bonds_rep(self):
        """ Create the geometry for the atomic structure bonds representation. 

            The bond representaion visualizes DNA atomic bonds as lines. The sugar, purine, and pyrimidine ring atoms 
            are displayed as solid polygons. A list of which points belong to which base is created so that the base 
            number of a picked point on the bond geometry can be determined. 
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("Create atomic structure bonds rep.")
        self._logger.debug("Number of residues %d " % len(self.molecule.residues))
        num_bond_points = 0
        bond_indexes = []
        bond_points = []
        num_plane_points = 0
        plane_counts = []
        plane_points = []
        plane_indexes = []
        num_res = len(self.molecule.residues)
        res_list = list(self.molecule.residues.values())
        for i in xrange(0,num_res): 
            residue_atoms = res_list[i]
            atom = residue_atoms["P"]
            self._get_residue_bond_geometry(bond_points, atom.res_name, residue_atoms)
            self._get_residue_plane_geometry(plane_counts, plane_points, atom.res_name, residue_atoms)
            if i != num_res-1:
                next_residue_atoms = res_list[i+1]
                self._add_residue_bond("P", residue_atoms, "O3'", next_residue_atoms, bond_points)
            # Add indexes into the bond and plane points arrays to map points to individual bases.
            num_bond_points = len(bond_points)
            bond_indexes.append(num_bond_points)
            num_plane_points = len(plane_counts)
            plane_indexes.append(num_plane_points)
            #self._logger.debug("Residue %d  number of bounds points %d " % (i,num_bond_points))
            #self._logger.debug("Residue %d  number of plane points %d " % (i,num_plane_points))
        #__for i in xrange(0,num_res)

        # Create the lines geometry for atom bonds.
        name = "BaseBonds:%s" %  self.strand_name
        geom = VisGeometryLines(name, bond_points)
        geom.color = self.color[:]
        geom.line_width = 2.0
        geom.entity_indexes = bond_indexes
        geom.selected_callback = self.select_bonds
        self.representations[VisAtomicStructureRepType.BONDS] = [geom]
        self.graphics.add_render_geometry(geom)
        # Create the polygon geometry for sugar and base planes.
        name = "BasePlanes:%s" %  self.strand_name
        geom = VisGeometryPolygon(name, plane_counts, plane_points)
        geom.color = self.color[:]
        geom.entity_indexes = plane_indexes
        geom.selected_callback = self.select_bonds
        self.representations[VisAtomicStructureRepType.BONDS].append(geom)
        self.graphics.add_render_geometry(geom)

    def select_bonds(self, geom, index):
        """ Process the selection of a bond representation. 

            Arguments:
                geom (VisGeometry): The geometry selected.
                index (int): The index into the geometry selected.
        """
        num_res = len(self.molecule.residues)
        res_list = list(self.molecule.residues.values())
        if (index >= 0) and (index < num_res):
            residue_atoms = res_list[index]
            atom = residue_atoms["P"]
            self._logger.info("Selected bonds %s  Residue sequence number %d  Base name %s" % (self.name, atom.res_seq_num, 
                atom.res_name ))
            self.print_info()
        #__if residue_num >= 0 and residue_num < num_res

    def create_check_rep(self):
        """ Create the geometry for the atomic structure check representation. 

            This visualization is used to display P atoms bond lengths that deviate from too much from an average length.
            The geometry will be a set of thick lines for only those bonds whose deviation is larger than a given tolerance. 
        """
        #self._logger.setLevel(logging.DEBUG)
        self._logger.debug("Create atomic structure check rep.")
        points = []
        s = self.scale
        for atom in self.molecule.atoms:
            if atom.element.strip() == 'P':
                point = [s*atom.coords[0], s*atom.coords[1], s*atom.coords[2]]
                points.append(point)
        #__for atom in self.molecule.atoms

        # Check P atom distances.
        verts = []
        entity_indexes = []
        lengths = []
        n = 2
        for i in xrange(0,len(points)-1):
            pt1 = points[i]
            pt2 = points[i+1]
            v = [pt1[j] - pt2[j] for j in xrange(0,3)]
            dist = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
            if abs(dist - VisAtomicStructure.average_p_bond_length) > VisAtomicStructure.p_bond_length_tol:
                verts.append(pt1)
                verts.append(pt2)
                lengths.append(dist)
                entity_indexes.append(n)
                n += 2
        #__for i in xrange(0,len(points)-1)

        # Create the geometry.
        show_vertices = True
        name = "Check:%s" %  self.strand_name
        geom = VisGeometryLines(name, verts)
        geom.line_width = 8.0
        geom.data.append(("lengths",lengths))    # add bond lengths data.
        geom.color = [0.8,0.8,0.4,1.0]
        geom.selected_callback = self.select_check   
        geom.entity_indexes = entity_indexes
        self.representations[VisAtomicStructureRepType.CHECK] = [geom]
        self.graphics.add_render_geometry(geom)

    def select_check(self, geom, index):
        """ Process the selection of a check representation. 

            Arguments:
                geom (VisGeometry): The geometry selected.
                index (int): The index into the geometry selected.
        """
        if index == None:
            return
        geom = self.representations[VisAtomicStructureRepType.CHECK][0]
        lengths = geom.data[0][1]
        self._logger.info("Selected check \'%s\'  P-P length %g " % (self.name, lengths[index]))
        #__if residue_num >= 0 and residue_num < num_res

    def _get_residue_bond_geometry(self, points, res_name, residue_atoms):
        """ Create the atomic bond geometry for the given residue.

            Arguments:
                points (List[float[3]]): The list of points to store the geometry.                
                res_name (String): The name of the residue.
                residue_atoms (List[Atom]): The list of atoms for the residue.
        """
        res_name = res_name.lower()
        # Add the backbone bonds.
        self._add_bonds(residue_atoms, VisDnaBonds.backbone, points)
        # Add the base bonds.
        if res_name == "da":
            self._add_bonds(residue_atoms, VisDnaBonds.ade, points)
        elif res_name == "dc":
            self._add_bonds(residue_atoms, VisDnaBonds.cyt, points)
        elif res_name == "dg":
            self._add_bonds(residue_atoms, VisDnaBonds.gua, points)
        elif res_name == "dt":
            self._add_bonds(residue_atoms, VisDnaBonds.thy, points)

    def _get_residue_plane_geometry(self, plane_counts, plane_points, res_name, residue_atoms):
        """ Create the atomic ring plane geometry for the given residue.

            Arguments:
                plane_counts (List[int]): The list of number of points for each polygon.
                plane_points (List[List[float]]): The list of points for each polygon.
                res_name (String): The residue name to create geometry for.
                residue_atoms (List[Atom]): The list of atoms for the residue.

            The geometry added are the coordinates of the polygons representing the surgar and base planes obtained
            from their ring atoms.
        """
        res_name = res_name.lower()
        # Add the sugar ring.
        self._add_planes(residue_atoms, VisDnaPlanes.sugar, plane_counts, plane_points)
        # Add the purine or pyrimidine ring(s).
        if res_name == "da":
            self._add_planes(residue_atoms, VisDnaPlanes.ade, plane_counts, plane_points)
        elif res_name == "dc":
            self._add_planes(residue_atoms, VisDnaPlanes.cyt, plane_counts, plane_points)
        elif res_name == "dg":
            self._add_planes(residue_atoms, VisDnaPlanes.gua, plane_counts, plane_points)
        elif res_name == "dt":
            self._add_planes(residue_atoms, VisDnaPlanes.thy, plane_counts, plane_points)

    def _add_planes(self, residue_atoms, planes, counts, points):
        """ Add the plane geometry for the given residue.

            Arguments:
                residue_atoms (List[Atom]): The list of atoms for the residue.
                planes (List[Tuples]): The list of atom names for a ring of atoms.
                counts (List[int]): The list of number of points for each polygon.
                points (List[List[float]]): The list of points for each polygon.
        """
        s = self.scale
        for plane in planes: 
            counts.append(len(plane))
            for atom_name in plane: 
                atom = residue_atoms[atom_name]
                points.append( [s*atom.coords[0], s*atom.coords[1], s*atom.coords[2]] )
        #__for plane in planes

    def _add_bonds(self, residue_atoms, bonds, points):
        """ Add the bond geometry for the given residue.

            Arguments:
                residue_atoms (List[Atom]): The list of atoms for the residue.
                bonds (List[Tuples]): The list of atom names for the residue bonds.
                points (List[List[float]]): The list of points for each bond (a line).
        """
        s = self.scale
        for bond in bonds: 
            if( bond[0] not in residue_atoms) or (bond[1] not in residue_atoms):
                continue
            atom1 = residue_atoms[bond[0]]
            atom2 = residue_atoms[bond[1]]
            points.append( [s*atom1.coords[0], s*atom1.coords[1], s*atom1.coords[2]] )
            points.append( [s*atom2.coords[0], s*atom2.coords[1], s*atom2.coords[2]] )

    def _add_residue_bond(self, atom_name1, residue_atoms1, atom_name2, residue_atoms2, points):
        """ Add the points representing a bond between atoms of two residues.

            Arguments:
                atom_name1 (String): The name of the atom in the 1st residue. 
                residue_atoms1 (String): The atoms of the 1st residue. 
                atom_name2 (String): The name of the atom in the 2nd residue. 
                residue_atoms2 (String): The atoms of the 2nd residue. 
                points (List[List[float]]): The list of points for each bond.
        """
        s = self.scale
        atom1 = residue_atoms1[atom_name1]
        atom2 = residue_atoms2[atom_name2]
        points.append( [s*atom1.coords[0], s*atom1.coords[1], s*atom1.coords[2]] )
        points.append( [s*atom2.coords[0], s*atom2.coords[1], s*atom2.coords[2]] )

class VisDnaBonds:
    """ This class defines the atom names for atomic bonds. """
    backbone = [
        # Phosphate group.
        ("P", "OP1"), ("P", "OP2"), ("P", "OP3"), ("P", "O5'"),
        # Sugar
        ("O5'", "C5'"), ("C5'", "C4'"), ("C4'", "C3'"), ("C3'", "C2'"), ("C2'", "C1'"), ("C1'", "O4'"),
        ("O4'", "C4'"), ("C3'", "O3'") ]

    # Bases.
    ade = [ ("C1'", "N9"), ("N9", "C8"), ("C8", "N7"), ("N7", "C5"), ("C5", "C4"), ("C4", "N9"), ("C5", "C6"),
            ("C6", "N1"), ("N1", "C2"), ("C2", "N3"), ("N3", "C4") ]

    gua = [ ("C1'", "N9"), ("N9", "C8"), ("C8", "N7"), ("N7", "C5"), ("C5", "C4"), ("C4", "N9"),
            ("C5", "C6"), ("C6", "N1"), ("N1", "C2"), ("C2", "N3"), ("N3", "C4"), ("C6", "O6"), ("C2", "N2") ]

    cyt = [ ("C1'", "N1"), ("N1", "C2"), ("C2", "N3"), ("N3", "C4"), ("C4", "C5"), ("C5", "C6"), ("C6", "N1"),
            ("C4", "N4"), ("C2", "O2") ]

    thy = [ ("C1'", "N1"), ("N1", "C2"), ("C2", "N3"), ("N3", "C4"), ("C4", "C5"), ("C5", "C6"), ("C6", "N1"),
            ("C4", "O4"), ("C5", "C7") ]

class VisDnaPlanes:
    """ This class defines the atom names for ring atoms. """
    sugar = [ ("C4'", "C3'", "C2'", "C1'", "O4'") ]
    ade = [ ("N9", "C8", "N7", "C5", "C4"), ("C5", "C6", "N1", "C2", "N3", "C4") ]
    gua = [ ("N9", "C8", "N7", "C5", "C4"), ("C5", "C6", "N1", "C2", "N3", "C4") ]
    cyt = [ ("N1", "C6", "C5", "C4", "N3", "C2") ]
    thy = [ ("N1", "C6", "C5", "C4", "N3", "C2") ]

