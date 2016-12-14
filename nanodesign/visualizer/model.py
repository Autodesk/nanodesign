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

"""This module is used to visualize a DNA structure derived from a DNA design.

   A DNA structure contains information about the virtual helices, strands,
   domains, base connectivity etc. derived from a DNA design. The VisModel class
   in this module manages the visualization of these entities using various
   visualization representations (i.e. the geometry used to display an entity).
"""
import os
import numpy as np
from .atomic_struct import VisAtomicStructure
from .cmd import VisCommand
from .extent import VisExtent
from .geometry import VisGeometryBox, VisGeometryCircle, VisGeometryPath, VisGeometryNumber, \
    VisGeometryPolygon, VisGeometryLines, VisGeometryCylinder 
from .graphics import VisGraphics 
from .helix import VisHelix
from .menu import VisMenu,VisMenuItem
from .strand import VisStrand

try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *

except ImportError as e:
    print "Could not import PyOpenGL."
    raise e

class VisModelRepType:
    """ This class defines the model visualization representation types. """
    UNKNOWN = 'unknown'
    BOUNDING_BOX  = 'Bounding box'
    GEOMETRY = 'Geometry'
    HELIX_NUMBERS = 'Virtual helix numbers'
    HELIX_PROJECTION = 'Virtual helix projection'

class VisModel(object):
    """ This class is used to visualize the entities (virtual helices, strands, atomistic model) of a DNA design structure. 

        Attributes:
            atomic_structure (AtomicStructure): The atomic structure of a DNA structure derived from a DNA design.
            atomic_structures (Dict[VisAtomicStructure]): The list of objects for atomic structure representations.
            atomic_structure_names (List[String]): The list of atomic structure names.
            cmd_file_name (String): The name of the input visualization commands file. This may be None.
            command (VisCommand): The command processing object. 
            commands (String): The string containing visualization commands from the command line. 
            dna_structure (DnaStructure): The DNA structure derived from a DNA design.
            extent (VisExtent): The extent of the DNA structure. 
            file_name (String): The design file name.
            graphics (VisGraphics): The visualization graphics object that manages the display of geometry for
                various entity representations.
            helices (Dict[VisHelix]): The list of objects for virtual helix representations. 
            helix_names (List[String]): The list of virtual helix names.
            menu (VisMenu): The menu object that manages the popup menu.
            name (String): The model name. This is the base name of the input file name. 
            strands (Dict[VisStrand]): The list of objects for strand representations.
            strand_names (List[String]): The list of strand names.
    """

    def __init__(self, file_name, cmd_file_name, commands, dna_structure, atomic_structure):
        """ Initialize a VisModel object.

            Arguments:
                file_name (String): The design file name.
                cmd_file_name (String): The name of the input visualization commands file. This may be None.
                commands (String): The string containing visualization commands from the command line. 
                dna_structure (DnaStructure): The DNA structure derived from a DNA design.
                atomic_structure (AtomicStructure): The atomic structure of a DNA structure derived from a DNA design.
        """
        self.name = os.path.basename(file_name)
        self.file_name = file_name
        self.cmd_file_name = cmd_file_name
        self.commands = commands
        self.dna_structure = dna_structure
        self.extent = VisExtent()
        self.command = VisCommand(self, cmd_file_name, commands)
        self.graphics = VisGraphics(self.name, self.command)
        self.menu = None 
        self.helices = {} 
        self.helix_names = [] 
        self.strands = {} 
        self.strand_names = [] 
        self.atomic_structure = atomic_structure 
        self.atomic_structures = {} 
        self.atomic_structure_names = [] 
        self.bounding_box_geometry = []
        self.helix_numbers_geometry = []
        self.helix_projection_geometry = []
        self.structure_geometry = []
        self.domains_temperature_range = None
        self._logger = logging.getLogger(__name__)
        self._set_extent()
        # Generate auxiliary data (e.g. domains) needed for certain visualizations.
        self.dna_structure.compute_aux_data()
        # Create helix, strand and atomic structure objects for visualization.
        self._create_helices()
        self._create_strands()
        self._create_atomic_structures()

    def _create_helices(self):
        """ Create a VisHelix object for each DnaHelix object. """
        dna_structure = self.dna_structure 
        helix_list = sorted(list(dna_structure.structure_helices_map.values()), key=lambda x: x.lattice_num)
        self._logger.info("Number of virtual helices %d" % (len(helix_list))) 
        self.helix_names.append(VisMenuItem.ALL)
        self.helix_names.append(VisMenuItem.NONE)
        for helix in helix_list:
            vis_helix = VisHelix(self, self.graphics, dna_structure, helix)
            self.helices[vis_helix.name] = vis_helix 
            self.helix_names.append(vis_helix.name)
    #_create_helices

    def _create_strands(self):
        """ Create a VisStrand object for each DnaStrand object. """
        dna_structure = self.dna_structure 
        # Create a  list of strands sorted by start helix and start position.
        strand_list = []
        for strand in dna_structure.strands:
            vis_strand = VisStrand(self, self.graphics, dna_structure, strand)
            self.strands[vis_strand.name] = vis_strand 
            strand_list.append(vis_strand)
        #__for strand in dna_structure.strands
        # Create a list of sorted strand names.
        strand_list.sort(VisStrand.compare)
        self._logger.info("Number of strands %d" % (len(strand_list))) 
        self.strand_names.append(VisMenuItem.ALL)
        self.strand_names.append(VisMenuItem.NONE)
        for strand in strand_list:
            self.strand_names.append(strand.name)
        #__for strand in strand_list
    #_create_strands

    def _create_atomic_structures(self):
        """ Create a VisAtomicStructure object for each atomic structure object. """
        if not self.atomic_structure:
            return
        dna_structure = self.dna_structure 
        self.atomic_structure_names.append(VisMenuItem.ALL)
        self.atomic_structure_names.append(VisMenuItem.NONE)
        molecules = self.atomic_structure.generate_structure_ss()
        #molecules = self.atomic_structure.generate_structure()
        self._logger.info("Number of molecules %d" % (len(molecules))) 
        id = 1
        atomic_struct_list = []
        for molecule in molecules: 
            #self._logger.info("Molecule %d  nchains %d  chainIDs %s " % (id, len(molecule.chains), str(list(molecule.chains))))
            atomic_struct = VisAtomicStructure(id, self, molecule, self.graphics) 
            self.atomic_structures[atomic_struct.strand_name] = atomic_struct
            atomic_struct_list.append(atomic_struct)
            id += 1
        #__for molecule in molecules
        # Create a list of sorted atomic structure names.
        atomic_struct_list.sort(VisAtomicStructure.compare)
        for atomic_struct in atomic_struct_list:
            self.atomic_structure_names.append(atomic_struct.strand_name)
    #__create_atomic_structures

    def start_interactive(self):
        """ Start the interactive visualization of the model. """
        self.graphics.set_extent(self.extent)
        cell_width = 0.1
        #self.graphics.initial_xform.rotate_x = -90.0 
        self.graphics.initialize_graphics()
        # Create the popup menu.
        self._create_menu()
        # Execute commands from a file or the command line.
        self.command.execute_file_cmds()
        self.command.execute_cmds()
        # Show the bounding box and helix numbers.
        self.command.generate_model_cmd(VisModelRepType.BOUNDING_BOX, "true")
        self.command.generate_model_cmd(VisModelRepType.HELIX_NUMBERS, "true")
        self.graphics.start_interactive()

    def _set_extent(self):
        """ Set the model extent from the DNA structure. """
        for helix in self.dna_structure.structure_helices_map.itervalues():
            for coord in helix.helix_axis_coords:
                self.extent.update(coord[0], coord[1], coord[2])
        xmin,xmax,ymin,ymax,zmin,zmax = self.extent.get_bounds()
        self._logger.info("Extent  xmin %f  ymin %f  zmin %f" % (xmin, ymin, zmin))
        self._logger.info("        xmax %f  ymax %f  zmax %f" % (xmax, ymax, zmax))

    def _create_menu(self):
        """ Create the popup menu for visualizing helices, strands, etc. """
        self.menu = VisMenu(self.command, self.helix_names, self.strand_names, self.atomic_structure_names)
        #self._logger.info("Helix names %s " % str(self.helix_names))
        #self._logger.info("Strand names %s " % str(self.strand_names))
        #self._logger.info("Atomic structure names %s " % str(self.atomic_structure_names))
        self.graphics.menu = self.menu 

    def show_helix(self, name, rep, attributes):
        """ Show a helix with the given representation. 

            Arguments:
                name (String): The name of the helix. 
                rep (String): The representation name from VisHelixRepType.
                attributes (List[String,Object)]: The list of attributes for the helix. 
        """
        show = None
        color = None

        # Process attributes.
        for key,value in attributes: 
            if key == 'show':
                show = value
            elif key == 'color':
                color = value
        #__for name,value in attibutes

        if name not in self.helix_names:
            self._logger.error("Unknown helix named \'%s\' " % name) 
            return

        # Show or hide all helices.
        if name == VisMenuItem.ALL:
            self._logger.info("Show all ") 
            display = False
            for i,helix in enumerate(self.helices.values()):
                helix.show(rep,show,display)

        # Show a helix named 'name'.
        else:
            helix = self.helices[name]
            if color != None:
                helix.set_color(rep, color)
            if show != None:
                helix.show(rep,show)
        self.graphics.display()

    def show_atomic_struct(self, name, rep, show):
        """ Show an atomic structure with the given representation.

            Arguments:
                name (String): The name of the atomic structure. 
                rep (String): The representation name from VisAtomicStructureRepType.
                show (bool): If true then show the geometry for the representation, else hide it.
        """
        if name not in self.atomic_structure_names:
            self._logger.error("Unknown atomic structure named \'%s\' " % name)
            return
        if name == VisMenuItem.ALL:
            display = False
            for atom_struct in self.atomic_structures.values():
                atom_struct.show(rep,show,display)
            self.graphics.display()
        # Show a atom struct named 'name'.
        else:
            atomic_struct = self.atomic_structures[name]
            atomic_struct.show(rep,show)

    def show_strand(self, name, rep, attributes):
        """ Show a strand with the given representation. 

            Arguments:
                name (String): The name of the strand . 
                rep (String): The representation name from VisStrandRepType.
                attributes (List[String,Object)]: The list of attributes for the helix. 
        """
        # Process attributes.
        show = None
        color = None
        line_width = None
        for key,value in attributes:
            if key == 'show':
                show = value
            elif key == 'color':
                color = value
            elif key == 'line_width':
                line_width = value
        #__for name,value in attibutes

        if name not in self.strand_names:
            self._logger.error("Unknown strand named \'%s\' " % name)
            return

        # Show or hide all strands.
        if name == VisMenuItem.ALL:
            self._logger.info("Show all ")
            display = False
            for strand in self.strands.values():
                strand.show(rep,show,display)
            self.graphics.display()
        # Show a strand named 'name'.
        else:
            strand = self.strands[name]
            if color != None:
                strand.set_color(rep, color)
            if line_width != None:
                strand.set_line_width(rep, line_width)
            if show != None:
                strand.show(rep,show)
            self.graphics.display()

    def show_bounding_box(self, show):
        """ Show a box bounding the model extent. """
        self._logger.info("Show bounding box %s " % str(show))
        if show and len(self.bounding_box_geometry) == 0:
            self._create_bounding_box()
        for geom in self.bounding_box_geometry:
             geom.visible = show
        self.graphics.display()

    def show_structure_geometry(self, show):
        """ Show the structure geometry. """
        if show and len(self.structure_geometry) == 0:
            self._create_structure_geometry()
        for geom in self.structure_geometry:
             geom.visible = show
        self.graphics.display()

    def show_helix_numbers(self, show):
        """ Show the virtual helix numbers. """
        if show and len(self.helix_numbers_geometry) == 0:
            self._create_helix_numbers()
        for geom in self.helix_numbers_geometry:
             geom.visible = show
        self.graphics.display()

    def show_helix_projections(self, show):
        """ Show the virtual helix projections. """
        if show and len(self.helix_projection_geometry) == 0:
            self._create_helix_projections()
        for geom in self.helix_projection_geometry:
             geom.visible = show
        self.graphics.display()

    def _create_bounding_box(self):
        # Create bounding box.
        cx, cy, cz = self.extent.get_center()
        dx, dy, dz = self.extent.get_widths()
        s = 1.0
        box_extent = VisExtent()
        box_extent.xmin = cx - s*dx
        box_extent.xmax = cx + s*dx
        box_extent.ymin = cy - s*dy
        box_extent.ymax = cy + s*dy
        box_extent.zmin = cz - s*dz
        box_extent.zmax = cz + s*dz
        name = "BoundingBox" 
        bbox = VisGeometryBox(name, box_extent)
        bbox.color = [0.5,0.5,0.5,1.0]
        bbox.line_width = 2.0
        self.graphics.add_render_geometry(bbox)
        self.bounding_box_geometry.append(bbox)

        counts = [4]
        points = [ [box_extent.xmin, box_extent.ymin, box_extent.zmin], 
                   [box_extent.xmax, box_extent.ymin, box_extent.zmin], 
                   [box_extent.xmax, box_extent.ymax, box_extent.zmin], 
                   [box_extent.xmin, box_extent.ymax, box_extent.zmin] ]

        name = "BoundingBox:bottom" 
        reverse_normals = True
        geom = VisGeometryPolygon("bob", counts, points, reverse_normals)
        geom.color = [0.8,0.8,0.8,1.0]
        #self.graphics.add_render_geometry(geom)

    def _create_helix_numbers(self):
        """ Create the geometry for virtual helix ends projected onto box and numbered. """
        cx, cy, cz = self.extent.get_center()
        dx, dy, dz = self.extent.get_widths()
        s = 1.0
        ymin = cy - s*dy

        radius = 1.0
        color = [0.5,0.5,0.5,1.0]
        for helix in self.helices.values():
            point1 = helix.vhelix.end_coordinates[0].copy()
            point2 = helix.vhelix.end_coordinates[1].copy()
            direction = [point2[i] - point1[i] for i in range(3)]
            # Create a circle.
            point1[1] = ymin 
            name = "HelixEnd:" + str(helix.id)
            geom = VisGeometryCircle(name, radius, point1, direction)
            geom.line_width = 1.0
            geom.color = color 
            self.graphics.add_render_geometry(geom)
            self.helix_numbers_geometry.append(geom)
            # Create the helix number.
            width = 0.2*radius
            number = str(helix.id)
            point1 = helix.vhelix.end_coordinates[0].copy()
            point1[1] = ymin 
            name = "HelixEndNumber:" + str(helix.id)
            geom = VisGeometryNumber(name, number, width, point1, direction)
            geom.color = [0.0,0.0,0.0,1.0] 
            geom.line_width = 1.0 
            self.graphics.add_render_geometry(geom)
            self.helix_numbers_geometry.append(geom)

    def _create_helix_projections(self):
        """ Create the geometry for virtual helices projected along coordinate axes. """ 
        cx, cy, cz = self.extent.get_center()
        dx, dy, dz = self.extent.get_widths()
        s = 1.0
        xmin = cx - s*dx
        zmin = cz - s*dz
        color = [0.5, 0.7, 0.9, 1.0] 
        for helix in self.helices.values():
            # Create a helix axis in y-z plane at xmin.
            point1 = helix.vhelix.end_coordinates[0].copy()
            point2 = helix.vhelix.end_coordinates[1].copy()
            points = [point1, point2]
            point1[0] = xmin 
            point2[0] = xmin 
            name = "HelixYZAxis:" + str(helix.id)
            geom = VisGeometryPath(name, points)
            geom.color = color 
            self.graphics.add_render_geometry(geom)
            self.helix_projection_geometry.append(geom)
            # Create a helix axis in x-y plane at zmin.
            point1 = helix.vhelix.end_coordinates[0].copy()
            point2 = helix.vhelix.end_coordinates[1].copy()
            points = [point1, point2]
            point1[2] = zmin 
            point2[2] = zmin 
            name = "HelixXYAxis:" + str(helix.id)
            geom = VisGeometryPath(name, points)
            geom.color = color 
            self.graphics.add_render_geometry(geom)
            self.helix_projection_geometry.append(geom)
        #__for helix in self.helices.values()

    def _create_structure_geometry(self):
        """ Create the structure geometry. 

            The structure geometry shows dsDNA regions as solid cylinders.
        """
        dna_structure = self.dna_structure
        radius = 1.15
        num_sides = 16
        angle_offset = 45.0
        helix_list = list(self.helices.values())
        self._logger.info("Number of virtual helices %d" % (len(helix_list)))
        for helix in helix_list:
            # Get the indexes into helix_axis_coords of dsDNA regions.
            boundary_points = helix.get_boundaries()
            verts = []
            for i,points in enumerate(boundary_points): 
                pt1 = points[0]
                pt2 = points[1]
                #self._logger.info("Pt1 %s  pt2 %s" % (str(pt1), str(pt2))) 
                verts.append(pt1)
                verts.append(pt2)
                name = "ModelGeometry:%d_%d" % (helix.vhelix.lattice_num,i+1) 
                geom = VisGeometryCylinder(name, radius, pt1, pt2, num_sides)
                self.graphics.add_render_geometry(geom)
                self.structure_geometry.append(geom)
            #__for pos in boundary_pos 
            name = "ModelGeometry:%d" % (helix.vhelix.lattice_num) 
            arrows = False
            geom = VisGeometryLines(name, verts, arrows)
            geom.line_width = 2.0
            geom.color = [0.8,0.0,0.8,1]
            self.graphics.add_render_geometry(geom)
            self.structure_geometry.append(geom)

        #__for helix in helix_list

    def get_domains_temperature_range(self):
        """ Get the domains temperature range. """
        if self.domains_temperature_range == None:
           tmin = None
           tmax = None
           for domain in self.dna_structure.domain_list:
               temp = domain.melting_temperature()
               if temp == -500.0:
                   continue
               if not tmin:
                   tmin = temp
                   vmax = temp
               elif temp < tmin:
                   tmin = temp
               elif temp > tmax:
                   tmax = temp
           #__for domain in domain_list
           self._logger.info("Domain temperature range min %g  max %g" % (tmin, tmax))
           self.domains_temperature_range = (tmin,tmax)
        #__if self.domains_temperature_range == None

        return self.domains_temperature_range[0], self.domains_temperature_range[1]
    #__def get_domain_temperature_range

    #__create_structure_geometry(self)

