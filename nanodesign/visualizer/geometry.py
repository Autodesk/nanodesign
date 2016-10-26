#!/usr/bin/env python
""" This module is used to store and render differnt geometric objects.

    The classes in this module are used to store and render the geometry of the various representations for
    visualizing DNA design structures. The classes are derived from the VisGeometry base class.

    Geometry can be seleced using a 3D line defined by a graphics picking operation. The intersection
    of the 3D line with the geometry is calculated numerically using functions (e.g comp_line_line_intersect)
    defined here. The intersection point of the 3D geometry is stored together with the index of the selected
    geometry element. For example, a VisGeometryLines geometry represents a set of N lines. Picking on this
    geometry returns an index (between 0 and N-1) identifying the line element selected.

    Groups of elements (i.e. vertices, lines, polygons) of a geometry can often be associated with a single 
    entity from a DNA structure (e.g. a group of lines representing the atomic bonds a base). An array of 
    index counts is used to map entities to groups of elements. For example, suppose that a VisGeometryPath 
    geometry that displays a DNA strand as a set of N connected lines using N+1 vertices. Each vertex in the 
    path represnts a single base in the strand. The array of index counts for this geometry would look like:
    entity_indexes = [0, 1, 2, ..., N-1]. 

    Note that the OpenGL per-vertex operations (e.g. glColor, glNormal, glVertex) used here are now deprecated 
    and that for large models rendering will be slow.
"""
from abc import ABCMeta, abstractmethod, abstractproperty
import os
import sys
import random
import numpy as np
from .extent import VisExtent
from math import sqrt,cos,sin,pi,acos

try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *

except ImportError as e:
    print "Could not import PyOpenGL."
    raise e


class VisGeometry(object):
    """ This is the VisGeometry base class.

        Attributes:
            data (List): The geometry generic data. 
            color (List[float]): The default geometry RGBA color. This is a list of four floats.
            entity_indexes (List[int]): The list of entiy indexes for the geometry. 
            extent (VisExtent): The geometry 3D extent.
            id (int): The geometry ID; just the count of the current number of geometries. 
            intersect_index (int): The last intersection index into the geometry.
            intersect_point (List[Float]): The last intersection point for the geometry.
            line_width (Float): The geometry default rendering line width.
            name (String): The geometry name.
            num_vertices (int): The number of verices in the geometry. 
            selected (bool): If true then the geometry has been selected.
            selected_callback (Method): The function to call when the geometry is selected.
            selected_entity (int): The last entity selected.
            selected_vertex (int): The index into the vertices array of the closest vertex to the selected point.
            vertices ((NumPy Nx3 ndarray[float]): The geometry vertices.
            visible (Bool): If True then the geometry is visible.
            transparent (Bool): If True then the geometry is transparent.
    """
    __metaclass__ = ABCMeta
    num_geometries = 0
    vertices_size = 0
    connectivity_size = 0

    def __init__(self, name):
        self.name = name
        self.id = VisGeometry.num_geometries
        self.extent = VisExtent()
        self.vertices = None 
        self.num_vertices = 0 
        self.data = [] 
        self.color = [1.0, 1.0, 1.0, 1.0]
        self.highlight_color = [1.0,0.8,0.8,1.0]
        self.visible = True
        self.transparent = False
        self.line_width = 1.0
        self.intersect_index = None
        self.intersect_point = None
        self.entity_indexes = None
        self.selected = False
        self.selected_entity = None
        self.selected_callback = None
        self.selected_vertex = None
        VisGeometry.num_geometries += 1
    
    @abstractmethod
    def render():
        """ Render geometry. """ 
        raise NotImplementedError

    @abstractmethod
    def intersect_line(self, point1, point2):
        """ Intersect the geometry with a line. """ 
        raise NotImplementedError

    def select_entity(self):
        """ Select an entity using an index into a geometry's data (e.g. vertices) calculated by intersecting a 
            geometry with a pick line. 

            The self.entity_indexes list stores the index range for each entity in the geometry. Searching this list will
            will determine the entity referenced by the index.
        """
        if self.entity_indexes == None:
            return

        if (self.intersect_point != None) and (self.intersect_index != None): 
            for entity,index in enumerate(self.entity_indexes):
                if self.intersect_index < index:
                    self.selected_entity = entity
                    break
            #__for entity,index in enumerate(self.entity_indexes)
        #__if self.intersect_point

    def update_stats(self, size_conn, size_verts):
        VisGeometry.vertices_size += size_verts 
        VisGeometry.connectivity_size = size_conn 
        #print("Number of geometries %d" % VisGeometry.num_geometries) 
        #print("Size of connectivity %d" % VisGeometry.connectivity_size) 
        #print("Size of vertices %d" % VisGeometry.vertices_size) 

class VisGeometryBox(VisGeometry):
    """ This class is used to display the outline of a 3D box. """ 
    def __init__(self, name, extent):
        """ Initialize a VisGeometryBox object. 
            Arguments:
                name (String): The geometry name.
                extent (VisExtent): The box extent. 
        """
        VisGeometry.__init__(self, name)
        xmin,xmax,ymin,ymax,zmin,zmax = extent.get_bounds()
        self.xmin = xmin 
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax 
        self._create_geometry(xmin,xmax,ymin,ymax,zmin,zmax)

    def _create_geometry(self,xmin,xmax,ymin,ymax,zmin,zmax):
        """ Create the box geometry. """
        vertices = np.zeros((24,3), dtype=float);
        n = 0
        # Minimum Z face.
        z = zmin
        vertices[n,:] = [xmin, ymin, z]; n += 1
        vertices[n,:] = [xmax, ymin, z]; n += 1
        vertices[n,:] = [xmax, ymin, z]; n += 1
        vertices[n,:] = [xmax, ymax, z]; n += 1
        vertices[n,:] = [xmax, ymax, z]; n += 1
        vertices[n,:] = [xmin, ymax, z]; n += 1
        vertices[n,:] = [xmin, ymax, z]; n += 1
        vertices[n,:] = [xmin, ymin, z]; n += 1

        # Maximum Z face.
        z = zmax
        vertices[n,:] = [xmin, ymin, z]; n += 1
        vertices[n,:] = [xmax, ymin, z]; n += 1
        vertices[n,:] = [xmax, ymin, z]; n += 1
        vertices[n,:] = [xmax, ymax, z]; n += 1
        vertices[n,:] = [xmax, ymax, z]; n += 1
        vertices[n,:] = [xmin, ymax, z]; n += 1
        vertices[n,:] = [xmin, ymax, z]; n += 1
        vertices[n,:] = [xmin, ymin, z]; n += 1

        # Bottom Y face.
        vertices[n,:] = [xmin, ymin, zmin]; n += 1
        vertices[n,:] = [xmin, ymin, zmax]; n += 1
        vertices[n,:] = [xmax, ymin, zmin]; n += 1
        vertices[n,:] = [xmax, ymin, zmax]; n += 1

        # Top Y face.
        vertices[n,:] = [xmax, ymax, zmin]; n += 1
        vertices[n,:] = [xmax, ymax, zmax]; n += 1
        vertices[n,:] = [xmin, ymax, zmin]; n += 1
        vertices[n,:] = [xmin, ymax, zmax]; n += 1

        self.vertices = vertices 
        self.num_vertices = 24 
        self.update_stats(0, len(self.vertices))

    def render(self):
        """ Render the box geometry. """
        if not self.visible:
            return
        glDisable(GL_LIGHTING);
        glLineWidth(self.line_width)
        glColor4fv(self.color)
        glBegin(GL_LINES)
        for i in xrange(0,self.num_vertices):
            glVertex3dv(self.vertices[i])
        glEnd()
        if self.intersect_point: 
            draw_cross(self.intersect_point)
        glEnable(GL_LIGHTING);

    def intersect_line(self, point1, point2):
        """ Intersect the geometry with a line. """ 
        self.intersect_index = None 
        self.intersect_point = None
        intersect_points = [] 
        intersect_indexes = [] 
        for i in xrange(0,self.num_vertices):
            line2_p1 = self.vertices[i]
            if i == self.num_vertices-1:
                line2_p2 = self.vertices[0]
            else:
                line2_p2 = self.vertices[i+1]
            ipt = comp_line_line_intersect(point1, point2, line2_p1, line2_p2)
            if ipt: 
                intersect_indexes.append(i)
                intersect_points.append(ipt)
        #__for i in xrange(0,self.num_vertices)
        if len(intersect_points) != 0:
            i,pt = get_closest_point(intersect_points, point1)
            self.intersect_point = pt
            self.intersect_index = intersect_indexes[i]
        return (self.intersect_point != None)

class VisGeometryLines(VisGeometry):
    """ This class is used to display pairs of lines. 

        Each pair of the geometry's vertices is used to display a line. Arrowheads may be drawn at the
        end of each line. 
    """
    def __init__(self, name, points, arrows=False):
        """ Initialize a VisGeometryLines object. 
            Arguments:
                name (String): The geometry name.
                points (List(List[Float]): The list of 3D points defining the lines endpoints. 
                arrows (bool): If true then create arrowheads at the end of each line. 
        """
        VisGeometry.__init__(self, name)
        self.num_arrow_vertices = 0
        self.arrow_vertices = None 
        self._create_geometry(points, arrows)

    def _create_geometry(self, points, arrows):
        """ Create the lines geometry. """
        self.num_vertices = len(points)
        self.vertices = np.zeros((self.num_vertices,3), dtype=float);
        for i,point in enumerate(points):
            self.vertices[i,:] = point

        # Create the arrowhead geometry.
        if arrows:
            num_arrow_pts = 4
            self.num_arrow_vertices = len(points) * 2 * num_arrow_pts 
            self.arrow_vertices = np.zeros((self.num_arrow_vertices,3), dtype=float)
            n = 0
            for i in xrange(0,len(points)/2):
                pt1 = points[2*i]
                pt2 = points[2*i+1]
                arrow_verts = generate_arrow_geometry(pt1, pt2, num_arrow_pts)
                for vert in arrow_verts:
                    self.arrow_vertices[n,:] = vert 
                    n += 1
            #__for i in xrange(0,len(points)-1)
            self.update_stats(0, len(self.arrow_vertices))
        #__if arrows
        self.update_stats(0, len(self.vertices))

    def intersect_line(self, point1, point2):
        """ Intersect the geometry with a line. """
        if not self.visible:
            return False
        self.intersect_index = None
        intersect_indexes = []
        self.intersect_point = None
        intersect_points = []
        self.selected_entity = None
        for i in xrange(0,self.num_vertices-1):
            line2_p1 = self.vertices[i]
            line2_p2 = self.vertices[i+1]
            ipt = comp_line_line_intersect(point1, point2, line2_p1, line2_p2)
            if ipt:
                intersect_indexes.append(i)
                intersect_points.append(ipt)
        #__for i in xrange(0,self.num_vertices)

        # Set the selected entity.
        if len(intersect_points) != 0:
            i,pt = get_closest_point(intersect_points, point1)
            self.intersect_point = pt
            self.intersect_index = intersect_indexes[i]
        self.select_entity()
        return (self.intersect_point != None)

    def render(self):
        """ Render the lines. """
        if not self.visible:
            return
        glDisable(GL_LIGHTING);
        glLineWidth(self.line_width)
        glColor4fv(self.color)
        glBegin(GL_LINES)
        for i in xrange(0,self.num_vertices):
            glVertex3dv(self.vertices[i])
        glEnd()

        # Render the lines arrowheads. 
        if self.num_arrow_vertices: 
            glBegin(GL_LINES)
            for i in xrange(0,self.num_arrow_vertices):
                glVertex3dv(self.arrow_vertices[i])
            glEnd()

        if self.intersect_point: 
            draw_cross(self.intersect_point)

        # Highlight the selected line.
        if self.selected:
            glColor4fv([1.0,1.0,1.0,1.0])
            glLineWidth(6.0)
            glBegin(GL_LINES)
            if self.selected_entity == 0:
                i1 = 0
            else:
                i1 = self.entity_indexes[self.selected_entity-1]
            i2 = self.entity_indexes[self.selected_entity]
            for i in xrange(i1,i2):
                glVertex3dv(self.vertices[i])
            glEnd()
        glEnable(GL_LIGHTING);

class VisGeometryPath(VisGeometry):
    """ This class is used to display a contiguous set of points as a solid path. 

        Attributes:
            select_vertex (bool): If true then select path vertices rather than path lines. 
            sphere_radius (Float): The radius of path vertex spheres.
            start_marker (bool): If true then display the start of the path with a sphere.
            start_marker_radius (Float): The radius of the sphere at the path start.
            start_sphere (VisGeometrySphere): The geometry object for the start sphere.
            vertex_spheres (List[VisGeometrySphere]): The list of geometry objects for vertex spheres. 
            sphere_selected (bool): If true then a sphere has been selected.
            bend_points (List[List[Float]]): The list of points where the path bends to cross over to an adjacent helix.
    """
    def __init__(self, name, points, show_vertices=False, show_arrows=False):
        """ Initialize a VisGeometryLines object. 

            Arguments:
                name (String): The geometry name.
                points (List(List[Float]): The list of 3D points defining the path. 
                show_vertices(bool): If true then display spheres at the path verticess. 
                show_arrows (bool): If true then display arrows at the mid-sections of bends in the path. 
        """
        VisGeometry.__init__(self, name)
        self.start_marker = False
        self.select_vertex = False
        self.start_marker_radius = 0.1
        self.start_sphere = None
        self.sphere_radius = 0.05
        self.vertex_spheres = []
        self.sphere_selected = False
        self.bend_points = []
        self._create_geometry(points, show_vertices, show_arrows)

    def _create_geometry(self, points, show_vertices, show_arrows):
        """ Create the geometry for the path. 

            Arguments:
                points (List(List[Float]): The list of 3D points defining the path. 
                show_vertices(bool): If true then display spheres at the path verticess. 
                show_arrows (bool): If true then display arrows at the mid-sections of bends in the path. 
        """
        self.num_vertices = len(points)
        self.vertices = np.zeros((self.num_vertices,3), dtype=float);
        for i,point in enumerate(points):
            self.vertices[i,:] = point
            if show_vertices:
                name = self.name + ":" + str(i)
                sphere = VisGeometrySphere(name, self.vertices[i], self.sphere_radius)
                sphere.color = self.color 
                self.vertex_spheres.append(sphere)
        #__for i,point in enumerate(points)

        # If showing arrows then calculate the bend points and
        # the vertices for arrow heads.
        if show_arrows:
            # If the path is not circular then set an offset to allow 
            # creating an arrow head at the end of the path.
            v = [points[0][i] - points[-1][i] for i in range(3)]
            if vector_mag(v) == 0.0:
                offset = -1
            else:
                offset = 0
            # Calculate bend points.
            self.bend_points.append(points[0])
            for i in xrange(0,len(points)-1):
                point1 = points[i]
                point2 = points[i+1]
                v = [point2[i] - point1[i] for i in range(3)]
                dist = vector_mag(v)
                if dist > 0.4:
                    self.bend_points.append(point1[:])
                    self.bend_points.append(point2[:])
            #__for i in xrange(0,len(points)-1)
            last_point = points[-1]
            self.bend_points.append(last_point[:])
            # Create arrow geometry. Arrows are placed at the midpoint between bend points. 
            self.num_arrow_vertices = len(self.bend_points) * 8 
            self.arrow_vertices = np.zeros((self.num_arrow_vertices,3), dtype=float)
            n = 0
            for i in xrange(0,len(self.bend_points)+offset):
                # Add an arrowhead at the end of a non-circular path. 
                if i == len(self.bend_points)-1:
                    cpt = self.bend_points[i]
                # Calculate the midpoint between bend points. 
                else:
                    point1 = self.bend_points[i]
                    point2 = self.bend_points[i+1]
                    dir = [point2[i] - point1[i] for i in range(3)]
                    dist = vector_mag(dir)
                    if dist == 0.0:
                        continue
                    u = vector_norm(dir)
                    v,w = compute_basis(u)
                    cpt = [point1[i] + 0.5*dist*u[i] for i in range(3)]
                #__if i == len(self.bend_points)-1
                scale = 0.2
                hs = 0.05*scale
                asc = 0.04
                # Create an arrow head.
                for j in xrange(0,3):
                    self.arrow_vertices[n,j]   = cpt[j] + scale*u[j]
                    self.arrow_vertices[n+1,j] = cpt[j] + hs*u[j] + asc*v[j]
                    self.arrow_vertices[n+2,j] = cpt[j] + scale*u[j]
                    self.arrow_vertices[n+3,j] = cpt[j] + hs*u[j] - asc*v[j]
                    self.arrow_vertices[n+4,j] = cpt[j] + scale*u[j]
                    self.arrow_vertices[n+5,j] = cpt[j] + hs*u[j] + asc*w[j]
                    self.arrow_vertices[n+6,j] = cpt[j] + scale*u[j]
                    self.arrow_vertices[n+7,j] = cpt[j] + hs*u[j] - asc*w[j]
                #__for j in xrange(0,3)
                n += 8
            #__for i in xrange(0,len(self.bend_points)/2)

        #self.update_stats(0, len(self.arrow_vertices))
        #self.update_stats(0, len(self.vertices))

    def intersect_line(self, point1, point2):
        """ Intersect the geometry with a line. """ 

        if not self.visible:
            return False

        self.intersect_index = None
        self.intersect_point = None
        intersect_indexes = []
        intersect_points = []

        # Intersect lines.
        for i in xrange(0,self.num_vertices-1):
            line2_p1 = self.vertices[i]
            line2_p2 = self.vertices[i+1]
            ipt = comp_line_line_intersect(point1, point2, line2_p1, line2_p2)
            if ipt:
                intersect_indexes.append(i)
                intersect_points.append(ipt)
        #__for i in xrange(0,self.num_vertices)

        # Intersect vertex spheres.
        num_ipt = len(intersect_points)
        self.sphere_selected = False
        if self.vertex_spheres:
            for i,sphere in enumerate(self.vertex_spheres):
                if sphere.intersect_line(point1, point2):
                     intersect_indexes.append(i)
                     intersect_points.append(sphere.center)
            #__for i,sphere in enumerate(self.vertex_spheres)
        #__if self.vertex_spheres

        # Set the selected entity.
        if len(intersect_points) != 0:
            i,pt = get_closest_point(intersect_points, point1)
            self.intersect_point = pt
            self.intersect_index = intersect_indexes[i]
            self.selected_entity = self.intersect_index 
            if i > num_ipt-1:
                self.sphere_selected = True
        #__if len(intersect_points) != 0

        # Determine the closest vertex to the selected point. We will use this
        # to select the nearest vertex to the pick point if needed.
        if len(intersect_points) != 0:
            if (self.selected_entity == 0) or (self.selected_entity == self.num_vertices-1):
                self.selected_vertex = self.selected_entity 
            else:
                i1 = self.entity_indexes[self.selected_entity]
                i2 = self.entity_indexes[self.selected_entity+1]
                pt1 = self.vertices[i1]
                pt2 = self.vertices[i2]
                v = [pt1[i] - self.intersect_point[i] for i in range(3)]
                d1 = vector_mag(v)
                v = [pt2[i] - self.intersect_point[i] for i in range(3)]
                d2 = vector_mag(v)
                if (d2 < d1):
                    self.selected_vertex = i2
                else:
                    self.selected_vertex = i1
        #__if len(intersect_points) != 0
        return (self.intersect_point != None)

    def render(self):
        """ Render the path geometry.  """

        if not self.visible:
            return

        # Render the path lines. 
        glDisable(GL_LIGHTING);
        glLineWidth(self.line_width)
        glColor4fv(self.color)
        glBegin(GL_LINE_STRIP)
        for i in xrange(0,self.num_vertices):
            glVertex3dv(self.vertices[i])
        glEnd()
        glEnable(GL_LIGHTING);

        # Render vertex spheres.
        if self.vertex_spheres:
            # Determine if a sphere has been selected and should be highlighted.
            if not self.selected or (self.selected_entity == None):
                highlight = False
            elif self.sphere_selected:
                highlight = True
            elif self.select_vertex and (self.selected_entity != self.num_vertices-1):
                i1 = self.entity_indexes[self.selected_entity]
                i2 = self.entity_indexes[self.selected_entity+1]
                pt1 = self.vertices[i1]
                pt2 = self.vertices[i2]
                v = [pt1[i] - self.intersect_point[i] for i in range(3)]
                d1 = vector_mag(v)
                v = [pt2[i] - self.intersect_point[i] for i in range(3)]
                d2 = vector_mag(v)
                if d2 < d1: 
                    self.selected_entity = i2
                highlight = True
                self.sphere_selected = True
            else:
                highlight = False

            for i,sphere in enumerate(self.vertex_spheres):
                if highlight and (i == self.selected_entity): 
                    continue
                sphere.color = self.color 
                sphere.render()
            #__for i,sphere in enumerate(self.vertex_spheres)

            if highlight: 
                sphere = self.vertex_spheres[self.selected_entity]
                sphere.color = self.highlight_color
                sphere.render()
        #__if self.vertex_spheres

        # Render the start sphere.
        if self.start_marker:
            if not self.start_sphere:
                name = self.name + "_start"
                self.start_sphere = VisGeometrySphere(name, self.vertices[0], self.start_marker_radius)
                self.start_sphere.color = self.color
            if self.selected and self.select_vertex and (self.selected_entity == 0):
                self.start_sphere.color = self.highlight_color
            else:
                self.start_sphere.color = self.color
            self.start_sphere.render()

        # Show the intersecion point.
        if self.intersect_point: 
            draw_cross(self.intersect_point)

        # Highlight selected entity. This will either be a line between two 
        # vertices or a region around the selected vertix. 
        if self.selected and (self.entity_indexes != None) and (self.selected_entity != None) and (not self.sphere_selected): 
            glColor4fv(self.highlight_color)
            glLineWidth(8.0)
            glBegin(GL_LINES)
            i1 = self.entity_indexes[self.selected_entity]
            i2 = self.entity_indexes[self.selected_entity+1]
            # Highlight a region around the selected vertix. 
            if self.select_vertex:
                pt1 = self.vertices[i1]
                pt2 = self.vertices[i2]
                if self.selected_vertex == 0:
                    v = [pt2[i] - pt1[i] for i in range(3)]
                    mpt = [pt1[i] + 0.5*v[i] for i in range(3)]
                    glVertex3dv(pt1)
                    glVertex3dv(mpt)
                elif self.selected_vertex == self.num_vertices-1:
                    v = [pt1[i] - pt2[i] for i in range(3)]
                    mpt = [pt2[i] + 0.5*v[i] for i in range(3)]
                    glVertex3dv(mpt)
                    glVertex3dv(pt2)
                else:
                    pt1 = self.vertices[self.selected_vertex-1]
                    pt2 = self.vertices[self.selected_vertex]
                    pt3 = self.vertices[self.selected_vertex+1]
                    v1 = [pt1[i] - pt2[i] for i in range(3)]
                    v2 = [pt3[i] - pt2[i] for i in range(3)]
                    mpt1 = [pt2[i] + 0.5*v1[i] for i in range(3)]
                    mpt2 = [pt2[i] + 0.5*v2[i] for i in range(3)]
                    glVertex3dv(mpt1)
                    glVertex3dv(pt2)
                    glVertex3dv(pt2)
                    glVertex3dv(mpt2)
                #__if self.selected_vertex == 0

            # Highlight a line between vertices. 
            else:
                for i in xrange(i1,i2+1):
                    glVertex3dv(self.vertices[i])
            #__if self.select_vertex
            glEnd()
        #__if self.selected and (self.entity_indexes != None) and (self.selected_entity != None) and (not self.sphere_selected)

        # Display arrows at the mid points between bends and markers at bend points.
        if self.bend_points:
            glColor4fv(self.color)
            glPointSize(3.0)
            glBegin(GL_POINTS)
            for i,point in enumerate(self.bend_points):
                if i != 0:
                    glVertex3dv(point)
            glEnd()
            glLineWidth(2.0)
            glDisable(GL_LIGHTING);
            glBegin(GL_LINES)
            for i in xrange(0,self.num_arrow_vertices):
                glVertex3dv(self.arrow_vertices[i])
            glEnd()
            glEnable(GL_LIGHTING);

    #__def render(self)

class VisGeometrySphere(VisGeometry):
    """ This class is used to display a solid sphere. """ 
    def __init__(self, name, center, radius, num_sides=20):
        """ Initialize a VisGeometryLines object. 

            Arguments:
                name (String): The geometry name.
                center (List[Float]): The center of the sphere.                
                radius (Float): The sphere radius. 
                num_sides (int): The number of sides (polygons) used to represent the sphere.
        """
        VisGeometry.__init__(self, name)
        self.center = [center[0], center[1], center[2]] 
        self.radius = radius
        self.num_sides = num_sides 

    def intersect_line(self, point1, point2):
        """ Intersect the sphere with a line. """ 
        x0 = point1[0]; y0 = point1[1]; z0 = point1[2]
        v1 = [point1[i] - point2[i] for i in range(3)]
        v2 = [point1[i] - self.center[i] for i in range(3)]
        xd = v1[0]; yd = v1[1]; zd = v1[2];
        xc = v2[0]; yc = v2[1]; zc = v2[2];
        a = xd*xd + yd*yd + zd*zd
        b = 2.0 * (xd*xc + yd*yc + zd*zc)
        c = xc*xc + yc*yc + zc*zc - self.radius*self.radius
        d = b*b - 4.0*a*c
        if d < 0.0:
            intersect = False
        else:
            intersect = True
        return intersect 

    def render(self):
        """ Render the sphere. """
        if not self.visible:
            return
        glColor4fv(self.color)
        glPushMatrix()
        try:
            glTranslatef(*self.center)
            glutSolidSphere(self.radius, self.num_sides, self.num_sides)
        finally:
            glPopMatrix()
    #__def render(self)

class VisGeometryCylinder(VisGeometry):
    """ This class is used to display a cylinder. 

        Attributes:
            axis (List[Float]): The cylinder axis; defined by point2-point1. 
            length (Float): The length of the cylinder; defined by |point1-point2|.
            normals ((NumPy Nx3 ndarray[float]): The cylinder polygons vertex normals. 
            num_tri (int): The number of triangles used to represent the cylinder. 
            point1 (List[Float]): The first endpoint defining the cylinder axis. 
            point2 (List[Float]): The second endpoint defining the cylinder axis. 
            radius (Float): The cylinder radius.
            tri_conn ((NumPy 3*num_tri ndarray[int]): The cylinder polygons connectivity. 
            unit_axis (List[Float]): The normalized cylinder axis.

        The cylinder orientation (axis) and length is defined by two points.
    """ 
    def __init__(self, name, radius, point1, point2, num_sides=20):
        """ Initialize a VisGeometryCylinder object. 

            Arguments:
                name (String): The geometry name.
                radius (Float): The cylinder radius.
                point1 (List[Float]): The first endpoint defining the cylinder axis. 
                point2 (List[Float]): The second endpoint defining the cylinder axis. 
                num_sides (int): The number of sides (polygons) used to represent the cylinder.
        """
        VisGeometry.__init__(self, name)
        self.radius = radius
        self.ends_highlight_color = None
        self.point1 = point1
        self.point2 = point2
        v = [point2[i] - point1[i] for i in range(3)]
        self.length = vector_mag(v)
        self.axis = None
        self.unit_axis = None
        self.num_sides = num_sides
        self.normals = None 
        self.num_tri = 0 
        self.tri_conn = None
        self._generate_cyl()

    def intersect_line(self, point1, point2):
        """ Intersect the cylinder with a line. """ 
        if not self.visible:
            return False
        self.intersect_index = None
        self.intersect_point = None
        ipt = self.intersect_cyl_line(point1, point2)
        if ipt: 
            self.intersect_index = 0
            self.intersect_point = ipt
            # Set the selected entity. A cylinder has only one entity.
            self.selected_entity = 0
        return (ipt != None)
    #__def intersect_line(self, point1, point2)

    def intersect_cyl_line(self, point1, point2):
        """ Calculate the intersection point between a line and the cylinder. 

            Arguments:
                point1 (List[Float]): The first endpoint defining the line . 
                point2 (List[Float]): The second endpoint defining the line. 

            Returns the intersection point. 
        """
        ipt = None
        c = self.point1
        d = self.point2
        v = [point2[i] - point1[i] for i in range(3)]
        h = vector_norm(self.axis)
        o = point1
        oc = [o[i] - c[i] for i in range(3)]
        odh = vector_dot(o, h)
        vdh = vector_dot(v, h)
        cdh = vector_dot(c, h)

        # Check for intersection with caps.
        cpt_in = False
        tc1 = 1e6

        if vdh != 0.0:
            co = [c[i] - o[i] for i in range(3)]
            codh = vector_dot(co, h)
            tc1 = codh / vdh
            cpt1 = [o[i] + tc1*v[i] for i in range(3)]
            vec = [cpt1[i] - c[i] for i in range(3)]
            vec_length = vector_mag(vec)

            if vec_length <= self.radius:
                cpt_in = True
                tc = tc1
                cpt = cpt1

            do = [d[i] - o[i] for i in range(3)]
            dodh = vector_dot(do, h)
            tc2 = dodh / vdh
            cpt2 = [o[i] + tc2*v[i] for i in range(3)]
            vec = [cpt2[i] - d[i] for i in range(3)]
            vec_length = vector_mag(vec)

            if vec_length <= self.radius: 
                if cpt_in: 
                    if tc2 < tc1:
                        tc = tc2
                        cpt = cpt2
                else:
                    cpt_in = True
                    tc = tc2
                    cpt = cpt2
        #__if vdh != 0.0

        if cpt_in:
            ipt = cpt

        # Check intersection with sides. 
        spt_in = False
        f = [oc[i] - (odh - cdh)*h[i] for i in range(3)]
        g = [v[i] - vdh*h[i] for i in range(3)]
        c2 = 2*vector_dot(f,g)
        c1 = vector_dot(g,g)
        c3 = vector_dot(f,f) - self.radius*self.radius
        desc = c2*c2 - 4.0*c1*c3

        if desc < 0.0:
            return ipt

        t1 = (-c2 + sqrt(desc)) / (2*c1)
        t2 = (-c2 - sqrt(desc)) / (2*c1)

        if (t1 < 0.0) and (t2 < 0.0):
            return ipt

        k1 = odh + t1*vdh - cdh
        k2 = odh + t2*vdh - cdh

        if (k1 <= self.length) and (k1 >= 0.0):
            t = t1
            spt_in = True
        else:
            t1 = -1.0

        if (k2 <= self.length) and (k2 >= 0.0):
            t = t2
            spt_in = True
        else:
            t2 = -1.0

        if (t1 != -1.0) and (t2 != -1.0):
            if t1 < t2: 
                t = t1
            else:
                t = t2

        if spt_in:
            if cpt_in and (tc < t): 
                ipt = cpt
            else:
                ipt = [o[i] + t*v[i] for i in range(3)]

        return ipt
    #__def intersect_cyl_line(self, point1, point2)

    def render(self):
        """ Render the cylinder. """
        if not self.visible:
            return
        axis = self.axis
        num_sides = self.num_sides
        num_tri = self.num_tri 
        verts = self.vertices 
        conn = self.tri_conn
        norms = self.normals
        n = 3
        conn_count = 0
        glEnable(GL_LIGHTING);
        #glShadeModel(GL_SMOOTH);
        glColor4fv(self.color)
        # Render the side polygons.
        for i in xrange(0,num_tri-2*num_sides):
            glBegin(GL_POLYGON)
            for j in xrange(0,n):
                k = conn[conn_count]
                nx = norms[k][0]
                ny = norms[k][1]
                nz = norms[k][2]
                glNormal3f(nx, ny, nz)
                x = verts[k][0]
                y = verts[k][1]
                z = verts[k][2]
                glVertex3f(x,y,z)
                conn_count += 1
            glEnd()
        #__for i in xrange(0,num)__

        # Render the cap polygons.
        for i in xrange(0,2*num_sides):
            if ((i % 2) == 0):
               s = -1.0
            else:
               s = 1.0
            glBegin(GL_POLYGON)
            for j in xrange(0,n):
                k = conn[conn_count]
                nx = s*axis[0]
                ny = s*axis[1]
                nz = s*axis[2]
                glNormal3f(nx, ny, nz)
                x = verts[k][0]
                y = verts[k][1]
                z = verts[k][2]
                glVertex3f(x,y,z)
                conn_count += 1
            glEnd()
        #__for i in xrange(0,num)__
        glDisable(GL_LIGHTING);

        # If the cylinder is selected then hightlight its boundng circles at its ends and 
        # several lines end-to-end on its surface.
        if self.selected and self.intersect_point:
            draw_cross(self.intersect_point)
            conn_count = self.num_tri-2*self.num_sides
            verts = self.vertices 
            conn = self.tri_conn
            conn_count = 0
            glColor4fv(self.highlight_color)
            glLineWidth(2.0)
            glBegin(GL_LINE_LOOP)
            n = 0
            for i in xrange(0,num_sides):
                x = verts[n][0]
                y = verts[n][1]
                z = verts[n][2]
                glVertex3f(x,y,z)
                n += 1
            glEnd()
            glBegin(GL_LINE_LOOP)
            for i in xrange(0,num_sides):
                x = verts[n][0]
                y = verts[n][1]
                z = verts[n][2]
                glVertex3f(x,y,z)
                n += 1
            glEnd()
            n = 0
            num_segs = 4
            dn = num_sides / num_segs
            glBegin(GL_LINES)
            for i in xrange(0,num_segs):
                x = verts[n][0]
                y = verts[n][1]
                z = verts[n][2]
                glVertex3f(x,y,z)
                x = verts[n+num_sides][0]
                y = verts[n+num_sides][1]
                z = verts[n+num_sides][2]
                glVertex3f(x,y,z)
                n += dn
            glEnd()

        # Highlight the ends of the cylinder.
        if self.ends_highlight_color != None:
            verts = self.vertices
            glColor4fv(self.ends_highlight_color)
            glLineWidth(1.5)
            n = 0
            for i in xrange(0,2):
                glBegin(GL_LINE_LOOP)
                for j in xrange(0,num_sides):
                    x = verts[n][0]
                    y = verts[n][1]
                    z = verts[n][2]
                    glVertex3f(x,y,z)
                    n += 1
                #__for j in xrange(0,num_sides)
                glEnd()
            #__for i in xrange(0,2)
    #__render()

    def _generate_cyl(self):
        """ Generate the cylinder polygons. """
        origin = self.point1
        self.axis = [self.point2[i] - self.point1[i] for i in range(3)]
        self.unit_axis = vector_norm(self.axis)
        u = self.unit_axis
        v,w = compute_basis(u)
        n = self.num_sides
        dt = 2.0*pi / n
        t = 0.0
        radius = self.radius 
        num_verts = 2*n + 2
        verts = np.zeros((num_verts, 3), dtype=float)
        vnorms = np.zeros((num_verts, 3), dtype=float)

        # Cylinder sides. 
        for i in xrange(0,n):
            ct = radius*cos(t)
            st = radius*sin(t)
            for j in xrange(0, 3):
                verts[i][j]    = origin[j] + ct*v[j] + st*w[j] 
                verts[i+n][j]  = origin[j] + ct*v[j] + st*w[j] + self.axis[j]
                vnorms[i][j]   = ct*v[j] + st*w[j]
                vnorms[i+n][j] = ct*v[j] + st*w[j]
            #__for j in xrange(0, 3)__
            t += dt
        #__for i in xrange(0,n)__

        # Cylinder end caps.
        for j in xrange(0, 3):
            verts[2*n][j]    = origin[j] 
            verts[2*n+1][j]  = origin[j] + self.axis[j]
            vnorms[2*n][j]   = -u[j];
            vnorms[2*n+1][j] = u[j];
        #__for j in xrange(0, 3)

        # Create the side polygons connectivity.
        num_tri = 4*n
        conn = np.zeros((3*num_tri), dtype=int)
        num_tri = 0;
        for i in xrange(0, n):
            i1 = i
            i2 = i+1
            i3 = i+n
            if (i2 == n):
                i2 = 0
            conn[3*num_tri+0] = i1
            conn[3*num_tri+1] = i2
            conn[3*num_tri+2] = i3
            num_tri += 1
            j1 = i2
            j2 = i2+n
            j3 = i3
            conn[3*num_tri+0] = j1
            conn[3*num_tri+1] = j2
            conn[3*num_tri+2] = j3
            num_tri += 1
        #__for i in xrange(0, n)__

        # Create the end caps polygons connectivity.
        for i in xrange(0,n):
            i1 = i+1 
            i2 = i 
            i3 = 2*n
            if (i1 == n): i1 = 0
            conn[3*num_tri+0] = i1
            conn[3*num_tri+1] = i2
            conn[3*num_tri+2] = i3
            num_tri += 1
            j2 = i+1+n 
            j1 = i+n 
            j3 = 2*n+1
            if (j2 == 2*n): j2 = n
            conn[3*num_tri+0] = j1
            conn[3*num_tri+1] = j2
            conn[3*num_tri+2] = j3
            num_tri += 1
        #__for i in xrange(0,n)

        self.vertices = verts
        self.normals = vnorms
        self.num_tri = num_tri 
        self.tri_conn = conn
    #_generate_cyl(self)

class VisGeometryAxes(VisGeometry):
    """ This class is used to display axes as three pairs of lines.

        An axes is draw as a set of three lines originating from a common origin. This is used to display 
        coordinate frames. An arrowhead is drawn for the third axis to better show directionality 
        (e.g. helix polarity). 
    """
    def __init__(self, name, origins, directions, scale=1.0):
        """ Initialize a VisGeometryAxes object. 

            Arguments:
                name (String): The geometry name.
                origins (List[List[Float]]): The list of N axis origins. 
                directions ((NumPy 3x3xN ndarray[float]): The axis three directions.  ith axis1 = directions[:,0,i], 
                    ith axis2 = directions[:,1,i], ith axis3 = directions[:,2,i]
                scale (Float): The axis length scale.
        """
        VisGeometry.__init__(self, name)
        self.line_width = 1.0
        self.color = [1.0, 0.0, 1.0, 1.0]
        self._create_geometry(origins, directions, scale)

    def _create_geometry(self, origins, directions, scale):
        """ Create the axes geometry. """
        num_axes = len(origins)
        self.num_vp = 14
        self.num_vertices = self.num_vp*num_axes
        self.vertices = np.zeros((self.num_vertices,3), dtype=float)
        n = 0
        hs = 0.8*scale
        asc = 0.1*hs
        for i in xrange(0,num_axes):
            for j in xrange(0,3):
                self.vertices[n,j] = origins[i,j] 
                self.vertices[n+1,j] = origins[i,j] + scale*directions[j,0,i]
                self.vertices[n+2,j] = origins[i,j] 
                self.vertices[n+3,j] = origins[i,j] + scale*directions[j,1,i]
                self.vertices[n+4,j] = origins[i,j] 
                self.vertices[n+5,j] = origins[i,j] + scale*directions[j,2,i]
                # Arrow head.
                self.vertices[n+6,j] = origins[i,j] + scale*directions[j,2,i]
                self.vertices[n+7,j] = origins[i,j] + hs*directions[j,2,i] + asc*directions[j,1,i]
                self.vertices[n+8,j] = origins[i,j] + scale*directions[j,2,i]
                self.vertices[n+9,j] = origins[i,j] + hs*directions[j,2,i] - asc*directions[j,1,i]
                self.vertices[n+10,j] = origins[i,j] + scale*directions[j,2,i]
                self.vertices[n+11,j] = origins[i,j] + hs*directions[j,2,i] + asc*directions[j,0,i]
                self.vertices[n+12,j] = origins[i,j] + scale*directions[j,2,i]
                self.vertices[n+13,j] = origins[i,j] + hs*directions[j,2,i] - asc*directions[j,0,i]
            #__for j in xrange(0,num_axes)
            n += self.num_vp
        #__for i in xrange(0,num_axes)

    def intersect_line(self, point1, point2):
        """ Intersect the geometry with a line. """ 
        if not self.visible:
            return False
        self.intersect_index = None
        self.intersect_point = None
        self.selected_entity = None
        for i in xrange(0,self.num_vertices-1):
            line2_p1 = self.vertices[i]
            line2_p2 = self.vertices[i+1]
            ipt = comp_line_line_intersect(point1, point2, line2_p1, line2_p2)
            if ipt:
                self.selected_entity = i/self.num_vp
                self.intersect_index = i
                self.intersect_point = ipt
                return True
        #__for i in xrange(0,self.num_vertices)
        return False 

    def render(self):
        """ Render the axes geometry. """
        if not self.visible:
            return
        glLineWidth(self.line_width)
        glColor4fv(self.color)
        glDisable(GL_LIGHTING);
        glBegin(GL_LINES)
        for i in xrange(0,self.num_vertices):
            glVertex3dv(self.vertices[i])
        glEnd()
        glEnable(GL_LIGHTING);

        # Highlight the selected axis.
        if self.selected and (self.selected_entity != None):
            glColor4fv(self.highlight_color)
            glLineWidth(6.0)
            glBegin(GL_LINES)
            i1 = self.selected_entity*self.num_vp
            i2 = i1 + self.num_vp 
            for i in xrange(i1,i2):
                glVertex3dv(self.vertices[i])
            glEnd()
    #__def render(self)

class VisGeometryPolygon(VisGeometry):
    """ This class is used to display polygons. 

        Attributes:
            centers ((NumPy Nx3 ndarray[float]): The polygons face centers. 
            counts (List[int]): The list of the number of points for each polygon. 
            normals ((NumPy Nx3 ndarray[float]): The polygons face normals. 
    """
    def __init__(self, name, counts, points, reverse_normals=False):
        """ Initialize a VisGeometryPolygon object. 

            Arguments:
                name (String): The geometry name.
                points (List(List[Float]): The list of 3D points defining the lines endpoints. 
        """
        VisGeometry.__init__(self, name)
        self.num_polygons = len(counts)
        self.num_vertices = len(points)
        self.vertices = np.zeros((self.num_vertices, 3), dtype=float)
        self.centers = np.zeros((self.num_polygons, 3), dtype=float)
        self.normals = np.zeros((self.num_polygons, 3), dtype=float)
        self.counts = counts[:]

        # Calculate the center of each polygon. 
        pindex = 0
        for i,n in enumerate(counts):
            cx = 0.0
            cy = 0.0
            cz = 0.0
            for j in xrange(0,n):
                self.vertices[pindex,:] = points[pindex]
                cx += points[pindex][0]
                cy += points[pindex][1]
                cz += points[pindex][2]
                pindex += 1
            #__for point in points
            self.centers[i] = [ cx/n, cy/n, cz/n ]
        #__for i,n in enumerate(counts)

        # Compute normals.
        pindex = 0
        verts = self.vertices
        if reverse_normals:
            s = -1.0
        else:
            s = 1.0
        for i,n in enumerate(self.counts):
            nx = 0.0
            ny = 0.0
            nz = 0.0
            for j in xrange(0,n):
                if j == n-1:
                    k = 0
                else:
                    k = j+1
                nx += (verts[j+pindex][1] - verts[k+pindex][1]) * (verts[j+pindex][2] + verts[k+pindex][2])  
                ny += (verts[j+pindex][2] - verts[k+pindex][2]) * (verts[j+pindex][0] + verts[k+pindex][0])  
                nz += (verts[j+pindex][0] - verts[k+pindex][0]) * (verts[j+pindex][1] + verts[k+pindex][1])  
            #__for j in xrange(0,n)
            mag = sqrt(nx*nx + ny*ny + nz*nz)
            self.normals[i] = [ s*nx/mag, s*ny/mag, s*nz/mag ]
            pindex += n
        #__for i,n in enumerate(counts)

    def intersect_line(self, point1, point2):
        """ Intersect polygons with a line. """ 
        self.intersect_index = None
        self.intersect_point = None
        intersect_indexes = []
        intersect_points = [] 
        pindex = 0
        for i,n in enumerate(self.counts):
            normal = self.normals[i]
            for j in xrange(0,n):
                if j == n-1:
                    k = 0
                else:
                    k = j+1
                poly_pts = [self.vertices[j+pindex], self.vertices[k+pindex], self.centers[i] ]
                ipt = comp_poly_line_intersect(normal, poly_pts, point1, point2)
                if ipt: 
                    intersect_indexes.append(i)
                    intersect_points.append(ipt)
            #__for j in xrange(0,n)
            pindex += n
        #__for i,n in enumerate(counts)

        # Set the selected entity.
        if len(intersect_points) != 0:
            i,pt = get_closest_point(intersect_points, point1)
            self.intersect_point = pt
            self.intersect_index = intersect_indexes[i]
            self.select_entity()
        return (self.intersect_point != None)

    def render(self):
        """ Render the polygons. """
        if not self.visible:
            return

        glEnable(GL_LIGHTING);
        glDisable(GL_CULL_FACE);
        glLineWidth(self.line_width)
        glColor4fv(self.color)
        glBegin(GL_TRIANGLES);

        # If an entity in the geomertry has been seleced get 
        # the range of indexes in the counts[] array. 
        if self.selected_entity != None:
            if self.selected_entity == 0:
                i1 = 0
            else:
                i1 = self.entity_indexes[self.selected_entity-1]
            i2 = self.entity_indexes[self.selected_entity]
            highlight = True
        else:
            highlight = False
        #__if self.selected_entity != None

        # Render polygons storing selected vertices for later highlighting.
        highlight_vindex = []
        pindex = 0
        for i,n in enumerate(self.counts):
            glNormal3dv(self.normals[i])
            for j in xrange(0,n):
                if j == n-1:
                    k = 0
                else:
                    k = j+1
                if highlight and (i >= i1) and (i < i2):
                   highlight_vindex.append([j+pindex, k+pindex, i])
                else:
                    glVertex3dv(self.vertices[j+pindex])
                    glVertex3dv(self.vertices[k+pindex])
                    glVertex3dv(self.centers[i])
            #__for j in xrange(0,n)
            pindex += n
        #__for i,n in enumerate(counts)

        # Highlight the selected polygon.
        if highlight:
            glColor4fv(self.highlight_color)
            for ivert in highlight_vindex:
                glVertex3dv(self.vertices[ivert[0]])
                glVertex3dv(self.vertices[ivert[1]])
                glVertex3dv(self.centers[ivert[2]])
        #__if highlight
        glEnd()
        glEnable(GL_CULL_FACE);
        glDisable(GL_LIGHTING);

        # Draw intersection points.
        if self.selected and self.intersect_point:
            draw_cross(self.intersect_point)

#__class VisGeometryPolygon

class VisGeometrySymbols(VisGeometry):
    """ This class is used to display 3D symbols.

    """
    BASE_INSERT = "base_insert"
    BASE_DELETE = "base_delete"

    def __init__(self, name, symbol, origins, directions, scale=1.0):
        """ Initialize a VisGeometryLines object. 

            Arguments:
                name (String): The geometry name.
                origins (List[List[Float]]): The list of N axis origins. 
                directions ((NumPy 3x3xN ndarray[float]): The symbol three directions.  ith axis1 = directions[:,0,i], 
                    ith axis2 = directions[:,1,i], ith axis3 = directions[:,2,i]
                scale (Float): The symbol scale.
        """
        VisGeometry.__init__(self, name)
        self.line_width = 1.0
        self.color = [1.0, 0.0, 1.0, 1.0]
        self.symbol = symbol
        self._create_geometry(origins, directions, scale)

    def _create_geometry(self, origins, directions, scale):
        """ Create the axes geometry. """
        if self.symbol == VisGeometrySymbols.BASE_INSERT:
            self._create_insert_geometry(origins, directions, scale)
        elif self.symbol == VisGeometrySymbols.BASE_DELETE:
            self._create_delete_geometry(origins, directions, scale)

    def _create_delete_geometry(self, origins, directions, scale):
        """ Create the geometry for the delete symbol, a 3D cross. """
        num_symbols = len(origins)
        self.num_vp = 6
        self.num_vertices = self.num_vp*num_symbols
        self.vertices = np.zeros((self.num_vertices,3), dtype=float)
        n = 0
        for i in xrange(0,num_symbols):
            for j in xrange(0,3):
                self.vertices[n,j]   = origins[i,j] - scale*directions[j,0,i]
                self.vertices[n+1,j] = origins[i,j] + scale*directions[j,0,i]
                self.vertices[n+2,j] = origins[i,j] - scale*directions[j,1,i]
                self.vertices[n+3,j] = origins[i,j] + scale*directions[j,1,i]
                self.vertices[n+4,j] = origins[i,j] - scale*directions[j,2,i]
                self.vertices[n+5,j] = origins[i,j] + scale*directions[j,2,i]
            #__for j in xrange(0,3)
            n += self.num_vp
        #__for i in xrange(0,num_symbols)

    def _create_insert_geometry(self, origins, directions, scale):
        """ Create the insert geometry, a 3D triangle. """
        num_symbols = len(origins)
        self.num_vp = 6
        self.num_vertices = self.num_vp*num_symbols
        self.vertices = np.zeros((self.num_vertices,3), dtype=float)
        n = 0
        for i in xrange(0,num_symbols):
            v1 = directions[:,2,i]
            v2 = directions[:,0,i]
            for j in xrange(0,3):
                self.vertices[n,j] = origins[i,j]
                self.vertices[n+1,j] = origins[i,j] + scale*v1[j] + scale*v2[j]
                self.vertices[n+2,j] = origins[i,j] 
                self.vertices[n+3,j] = origins[i,j] - scale*v1[j] + scale*v2[j]
                self.vertices[n+4,j] = origins[i,j] - scale*v1[j] + scale*v2[j]
                self.vertices[n+5,j] = origins[i,j] + scale*v1[j] + scale*v2[j]
            #__for j in xrange(0,3)
            n += self.num_vp
        #__for i in xrange(0,num_symbols)

    def intersect_line(self, point1, point2):
        """ Intersect the geometry with a line. """ 
        if not self.visible:
            return False
        self.intersect_index = None
        self.intersect_point = None
        self.selected_entity = None
        for i in xrange(0,self.num_vertices-1):
            line2_p1 = self.vertices[i]
            line2_p2 = self.vertices[i+1]
            ipt = comp_line_line_intersect(point1, point2, line2_p1, line2_p2)
            if ipt:
                self.selected_entity = i/self.num_vp
                self.intersect_index = i
                self.intersect_point = ipt
                return True
        #__for i in xrange(0,self.num_vertices)
        return False 

    def render(self):
        """ Render the axes geometry. """
        if not self.visible:
            return
        glLineWidth(self.line_width)
        glColor4fv(self.color)
        glDisable(GL_LIGHTING);
        glBegin(GL_LINES)
        for i in xrange(0,self.num_vertices):
            glVertex3dv(self.vertices[i])
        glEnd()
        glEnable(GL_LIGHTING);

        # Highlight the selected axis.
        if self.selected and (self.selected_entity != None):
            glColor4fv(self.highlight_color)
            glLineWidth(6.0)
            glBegin(GL_LINES)
            i1 = self.selected_entity*self.num_vp
            i2 = i1 + self.num_vp 
            for i in xrange(i1,i2):
                glVertex3dv(self.vertices[i])
            glEnd()
    #__def render(self)

class VisGeometryCircle(VisGeometry):
    """ This class is used to display a circle in 3D. 

        Attributes:
            center (List[Float]): The circle center.                
            normal (List[Float]): The normal to the plane the circle lies in.                
            radius (Float): The circle radius. 
    """

    def __init__(self, name, radius, center, normal, num_sides=20):
        """ Initialize a VisGeometryLines object. 

            Arguments:
                name (String): The geometry name.
                radius (Float): The circle radius. 
                center (List[Float]): The circle center.                
                normal (List[Float]): The normal to the plane the circle lies in.                
                num_sides (int): The number of points used to represent the circle.
        """
        VisGeometry.__init__(self, name)
        self.radius = radius
        self.center = center
        self.normal = normal
        self.axis_normal = None
        self.num_sides = num_sides
        self._generate_circle()

    def intersect_line(self, point1, point2):
        """ Intersect the circle with a line. """ 
        ipt = None
        c = self.center
        v = [point2[i] - point1[i] for i in range(3)]
        h = vector_norm(self.normal)
        vdh = vector_dot(v, h)

        # Check for intersection when line and circle plane are not perpendicular.
        if vdh != 0.0:
            co = [c[i] - point1[i] for i in range(3)]
            codh = vector_dot(co, h)
            tc = codh / vdh
            cpt = [point1[i] + tc*v[i] for i in range(3)]
            vec = [cpt[i] - c[i] for i in range(3)]
            vec_length = vector_mag(vec)
            if vec_length <= self.radius:
                ipt = cpt 
        #__if vdh != 0.0

        # Check for intersection when line and circle plane are perpendicular.
        # Find the intersection by intersecting a line passing through the circle
        # center and perpendicular to the pick line. Then solve for the angle
        # for the intersected point to get the point on the circle.
        if ipt == None: 
            vn = vector_norm(v)
            u = cross_product(self.axis_normal,vn)
            cpt1 = [c[i] + self.radius*u[i] for i in range(3)]
            cpt2 = [c[i] - self.radius*u[i] for i in range(3)]
            ipt = comp_line_line_intersect(cpt1, cpt2, point1, point2)
            if ipt != None:
                vr = [c[i] - ipt[i] for i in range(3)]
                dp = vector_dot(vr, u)
                s = np.sign(dp)
                a = vector_mag(vr)
                angle = acos(a/self.radius)
                b = self.radius*sin(angle)
                ipt = [c[i] - s*a*u[i] + b*vn[i] for i in range(3)]
        #__if vdh != 0.0

        if ipt != None:
            self.intersect_index = 0
            self.intersect_point = ipt
            # Set the selected entity. A circle has only one entity.
            self.selected_entity = 0

        return (ipt != None)

    def render(self):
        """ Render the circle. """
        if not self.visible:
            return
        glLineWidth(self.line_width)
        if self.selected:
            glColor4fv(self.highlight_color)
        else:
            glColor4fv(self.color)
        glBegin(GL_LINE_STRIP)
        for i in xrange(0,self.num_vertices):
            glVertex3dv(self.vertices[i])
        glEnd()

    def _generate_circle(self):
        """ Generate the circle geometry. """
        origin = self.center
        self.axis_normal = vector_norm(self.normal)
        u = self.axis_normal
        v,w = compute_basis(u)
        n = self.num_sides
        dt = 2.0*pi / (n-1)
        t = 0.0
        radius = self.radius
        num_verts = n
        verts = np.zeros((num_verts, 3), dtype=float)
        for i in xrange(0,n):
            ct = radius*cos(t)
            st = radius*sin(t)
            for j in xrange(0, 3):
                verts[i][j] = origin[j] + ct*v[j] + st*w[j]
            #__for j in xrange(0, 3)__
            t += dt
        #__for i in xrange(0,n)__
        self.num_vertices = num_verts
        self.vertices = verts
    #__generate_circle(self):
#__class VisGeometryCircle(VisGeometry):

class VisGeometryNumber(VisGeometry):
    """ This class is used to display a number in 3D.  

        Attributes:
            center (List[Float]): The point in 3D to center the number.
            normal (List[Float]): The normal to the plane to display the number in.
            number (String): The number to display.
            u (List[Float]): The x-axis in the normal plane to display the number in.
            v (List[Float]): The y-axis in the normal plane to display the number in.
            width (Float): The width of the number to display.

        This class implements a simple number drawing function to replace the glutStrokeCharacter function 
        which is not supported on all platforms.
    """
    # Data describing the line segment coordinates on a 3x3 raster for digits 0-9.
    number_patterns = { 
        '0' : [ (0,0), (1,0), (1,2), (0,2), (0,0) ],
        '1' : [ (1,0), (1,2)],
        '2' : [ (0,2), (1,2), (1,0), (2,0)],
        '3' : [ (0,2), (1,2), (1,1), (1,0), (0,0), (), (0,1),(1,1)],
        '4' : [ (0,2), (0,1), (1,1), (1,0), (), (1,2), (1,0)],
        '5' : [ (1,2), (0,2), (0,1), (1,1), (1,0), (0,0)],
        '6' : [ (0,2), (0,0), (1,0), (1,1), (0,1)],
        '7' : [ (0,2), (1,2), (1,0)],
        '8' : [ (0,0), (1,0), (1,2), (0,2), (0,0), (), (0,1), (1,1)],
        '9' : [ (1,1), (0,1), (0,2), (1,2), (1,0)] 
    }

    def __init__(self, name, number, width, center, normal, u=[1,0,0], v=[0,0,1]):
        """ Initialize a VisGeometryLines object. 

            Arguments:
                name (String): The geometry name.
                number (String): The number to display.
                width (Float): The width of the number to display.
                center (List[Float]): The point in 3D to center the number.
                normal (List[Float]): The normal to the plane to display the number in.
                u (List[Float]): The x-axis in the normal plane to display the number in.
                v (List[Float]): The y-axis in the normal plane to display the number in.
        """
        VisGeometry.__init__(self, name)
        self.number = number 
        self.width = width
        self.center = center
        self.normal = vector_norm(normal)
        self.u = vector_norm(u)
        self.v = vector_norm(v)
        self._generate_number()

    def intersect_line(self, point1, point2):
        """ Intersect the number with a line. """ 
        for i in xrange(0,self.num_vertices-1):
            pt1 = self.vertices[i]
            pt2 = self.vertices[i+1]
            ipt = comp_line_line_intersect(pt1, pt2, point1, point2)
            if ipt:
                self.intersect_index = i
                self.intersect_point = ipt
        #__for i in xrange(0,self.num_vertices)
        return (self.intersect_point != None)

    def render(self):
        """ Render the number. """
        if not self.visible:
            return
        glLineWidth(self.line_width)
        if self.selected:
            glColor4fv(self.highlight_color)
        else:
            glColor4fv(self.color)
        glBegin(GL_LINES)
        for i in xrange(0,self.num_vertices):
            glVertex3dv(self.vertices[i])
        glEnd()

    def _generate_number(self):
        """ Generate the number geomatry. """
        origin = self.center
        w = self.normal
        u = self.u
        v = self.v
        width = self.width
        nd = len(self.number) / 2.0
        offset = [origin[j] - nd*width*u[j] for j in xrange(0,3)]
        points = []
        n = 0
        for digit in self.number:
            pattern = VisGeometryNumber.number_patterns[digit]
            for i in xrange(0,len(pattern)-1):
                disp1 = pattern[i]
                disp2 = pattern[i+1]
                if not disp1 or not disp2:
                    continue
                point1 = [offset[j] + width*disp1[0]*u[j] + width*disp1[1]*v[j] for j in xrange(0,3)]
                point2 = [offset[j] + width*disp2[0]*u[j] + width*disp2[1]*v[j] for j in xrange(0,3)]
                points.append(point1)
                points.append(point2)
            n += 1
            offset = [offset[j] + 2.2*n*width*u[j] for j in xrange(0,3)]
        #__for digit in self.number
        self.num_vertices = len(points)
        self.vertices = np.array(points, dtype=float)
    #__def _generate_number(self)
#__class VisGeometryNumber(VisGeometry)

def generate_arrow_geometry(point1, point2, num_arrow_pts=4):
    """ Generate a set of points representing an arrowhead at the end of a line. 

        Arguments:
            point1 (List[Float]): The start point defining of the line . 
            point2 (List[Float]): The end point of the line. 
            num_arrow_pts (int): The number of points used for the arrowhead.

        Returns:
            verts (List[Float]): The points representing an arrowhead. 

        The arrowdhead is defined by pairs of points connecting to point2.
    """
    u = [point2[i] - point1[i] for i in range(3)]
    v,w = compute_basis(u)
    dt = 2.0*pi / num_arrow_pts
    t = 0.0
    h = 0.9
    h = 0.7
    r = 0.05
    center = [point1[i] + h*u[i] for i in range(3)]
    verts = [] 
    for i in xrange(0,num_arrow_pts):
       ct = r*cos(t)
       st = r*sin(t)
       verts.append(point2)
       verts.append([center[j] + ct*v[j] + st*w[j] for j in range(3)])
       t += dt
    #__for i in xrange(0,n)__
    return verts 

def draw_cross(point,width=0.08,color=[1.0,0.8,0.8,1.0],line_width=2.0):
    """ Draw a 3D cross at the given point. """
    return
    glColor4fv(color)
    glLineWidth(line_width)
    glBegin(GL_LINES)
    pt1 = [ point[0]-width, point[1], point[2] ]
    pt2 = [ point[0]+width, point[1], point[2] ]
    glVertex3dv(pt1)
    glVertex3dv(pt2)
    pt1 = [ point[0], point[1]-width, point[2] ]
    pt2 = [ point[0], point[1]+width, point[2] ]
    glVertex3dv(pt1)
    glVertex3dv(pt2)
    pt1 = [ point[0], point[1], point[2]-width ]
    pt2 = [ point[0], point[1], point[2]+width ]
    glVertex3dv(pt1)
    glVertex3dv(pt2)
    glEnd()

def get_closest_point(points, point):
    """ Get the closest point in an array of points to the given point. """
    if len(points) == 1:
        return 0,points[0]
    min_dist = None
    min_pt = None
    min_i = None
    for i,pt in enumerate(points):
        v = [point[j] - pt[j] for j in xrange(0,3)]
        dist = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
        if (min_dist == None) or (dist < min_dist):
            min_dist = dist
            min_pt = pt
            min_i = i
        #__if (min_dist == None) or (dist < min_dist)
    #__for i,pt in enumerate(points)
    return min_i, min_pt

#----------------------------------------------------#
#----------------- vector functions -----------------#
#----------------------------------------------------#

def compute_basis(normal):
    """ Compute an orthonormal basis for a vector. """
    u = [0.0,0.0,0.0]
    v = [0.0,0.0,0.0]
    u[0] = -normal[1]
    u[1] =  normal[0]
    u[2] =  0.0

    if ((u[0] == 0.0) and (u[1] == 0.0)):
      u[0] = 1.0

    mag = vector_mag(u)
    if (mag == 0.0):
      return

    for i in xrange(0,3):
       u[i] = u[i] / mag

    v = cross_product(normal,u)
    mag = vector_mag(v)
    if (mag != 0.0):
       for i in xrange(0,3):
          v[i] = v[i] / mag
    return u,v

def cross_product(u, v):
    return [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]

def vector_mag(u):
    return sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2])

def vector_dot(u, v):
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def vector_norm(u):
    mag = vector_mag(u)
    if (mag != 0.0):
       return [u[0]/mag,u[1]/mag,u[2]/mag]
    else:
       return [0.0,0.0,0.0]

#------------------------------------------------------#
#----------------- geometry functions -----------------#
#------------------------------------------------------#

def comp_line_line_intersect(line1_p1, line1_p2, line2_p1, line2_p2):
    """ Compute the intersection of two lines. 

        Arguments:
            line1_p1 (List[Float]): The start point of the first line . 
            line1_p2 (List[Float]): The end point of the first line . 
            line2_p1 (List[Float]): The start point of the second line . 
            line2_p2 (List[Float]): The end point of the second line . 

        Returns:
            ipt (List[Float]): The intersection point. 
    """ 
    tol = 0.1
    v1 = [line1_p2[i] - line1_p1[i] for i in range(3)]
    v2 = [line2_p2[i] - line2_p1[i] for i in range(3)]
    a = vector_dot(v1, v1)
    b = -vector_dot(v1, v2)
    c = -b
    d = -vector_dot(v2, v2)
    e = sum([v1[i]*line2_p1[i] - v1[i]*line1_p1[i] for i in range(3)])
    f = sum([v2[i]*line2_p1[i] - v2[i]*line1_p1[i] for i in range(3)])
    det = a*d - c*b

    if det == 0.0: 
        return None

    t = (e*d - f*b) / det
    u = (a*f - e*c) / det

    if (t < 0.0) or (t > 1.0): 
        return None

    if (u < 0.0) or (u > 1.0): 
        return None

    ipt = [0.0]*3
    dist = 0.0
    for i in xrange(0,3):
        s1 = line1_p1[i] + t*v1[i]
        s2 = line2_p1[i] + u*v2[i]
        dist += (s1 - s2) * (s1 - s2)
        ipt[i] = s2

    dist = sqrt(dist)
    if dist > tol:
        ipt =  None

    return ipt

#__def comp_line_line_intersect

def comp_poly_line_intersect(normal, poly_pts, line_p1, line_p2):
    """ Compute the intersection of a line and a polygon.

        Arguments:
            normal (List[Float]): The normal to the polygon.
            poly_pts (List[List[Float]]): The polygon points.
            line_p1 (List[Float]): The start point of the line.
            line_p2 (List[Float]): The end point of the line.

        Returns:
            ipt (List[Float]): The intersection point.
    """
    ipt = [0.0]*3
    tol = 0.0001
    dir = [line_p1[i] - line_p2[i] for i in range(3)]
    ndp = vector_dot(dir, normal)

    if abs(ndp) < tol: 
      return None

    pt0 = poly_pts[0]
    v = [pt0[i] - line_p2[i] for i in range(3)]
    t = vector_dot(v,normal) / ndp

    if (t + tol < 0) or (t-tol > 1.0): 
        return None

    ipt[0] = line_p2[0] + t*dir[0]
    ipt[1] = line_p2[1] + t*dir[1]
    ipt[2] = line_p2[2] + t*dir[2]
    i0, i1, i2 = vector_max_proj(normal)

    if not point_in_triangle(poly_pts, ipt, i0, i1, i2):
        return None
    return ipt 

def vector_max_proj(normal):
    """ Determine a vector's maximum projection (i.e. the indices of a vector's maximum components) . 

        Arguments:
            normal (List[Float]): The normal to find the maximum projection. 

        Returns:
            i1, i2, i3 (int): The maximum projections. 
    """
    nx = normal[0]
    ny = normal[1]
    nz = normal[2]
    proj = [abs(nx), abs(ny), abs(nz)]
    max_proj = max(proj)
    axis = proj.index(max_proj)
    if axis == 0: 
        i0 = 0
        i1 = 1
        i2 = 2
    elif axis == 1: 
        i0 = 1
        i1 = 0
        i2 = 2
    elif axis == 2: 
        i0 = 2
        i1 = 0
        i2 = 1
    return i0, i1, i2

def point_in_triangle(tri, point, i0, i1, i2):
    """ Determine if a point is in a triangle. 

        Arguments:
            tri (List[Float]): The list of triangle points. 
            point (List[Float]): The point to test. 
            i1, i2, i3 (int): The maximum projections of the polygon's normal.

        Returns:
            in_tri (bool): If true then the point is in the triangle. 
    """
    u0 = point[i1] - tri[0][i1]
    v0 = point[i2] - tri[0][i2]
    u1 = tri[1][i1] - tri[0][i1]
    v1 = tri[1][i2] - tri[0][i2]
    u2 = tri[2][i1] - tri[0][i1]
    v2 = tri[2][i2] - tri[0][i2]
    in_tri = False

    if u1 == 0.0: 
        b = u0 / u2
        if (b >= 0.0) and (b <= 1.0) and (v1 != 0.0): 
            a = (v0 - b*v2) / v1
            in_tri = (a >= 0.0) and ((a + b) <= 1.00001)
        #__if (b >= 0.0) and (b <= 1.0) and (v1 != 0.0)
  
    else:
        b = (v0*u1 - u0*v1) / (v2*u1 - u2*v1)
        if (b >= 0.0) and (b <= 1.00001): 
            a = (u0 - b*u2) / u1
            in_tri = (a >= 0.0) and ((a + b) <= 1.00001)
        #__if (b >= 0.0) and (b <= 1.00001)
    #__if u1 == 0.0

    return in_tri


