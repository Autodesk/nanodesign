#!/usr/bin/env python
""" This module manages the OpenGL graphics used to display the geometry created by visualization.

    The points, lines and polygons of the geometry created by various visualization representations 
    of a DNA design are displayed in a graphics window using the PyOpenGL graphics API. The classes
    in this module are used to manage viewing, window, mouse and keyboard events and picking of displayed
    geometry. The VisGraphics class stores all of the geometry objects that are rendered in a scene. The
    actual drawing of graphics primitives is performed in the geometry.py module.

    Note that the OpenGL per-vertex operations (e.g. glColor, glNormal, glVertex) used here are now deprecated 
    and that for large models rendering will be slow. 
"""
import copy
from math import ceil, sqrt
import os
import sys
import random
from .menu import VisMenu 
from .extent import VisExtent

try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *

except ImportError as e:
    print "Could not import PyOpenGL."
    raise e

class MouseActions:
   """ This class defines mouse action types. """
   NONE         = -1
   ROTATE_XY    = 0
   ROTATE_Z     = 1
   TRANSLATE_XY = 2
   TRANSLATE_Z  = 3
   ZOOM         = 4
   PICK         = 5

class VisGraphicsXform(object):
    """ This class stores the transformation for the graphics scene. 

        Attributes:
            translate_x (Float): The translation in the x direction.
            translate_y (Float): The translation in the y direction.
            translate_z (Float): The translation in the z direction.
            rotate_x (Float): The rotation about the x axis. Units it degrees.
            rotate_y (Float): The rotation about the y axis. Units it degrees.
            rotate_z (Float): The rotation about the z axis. Units it degrees.
            scale (Float): The zoom scaling factor.
    """
    def __init__(self):
        """ Initialize a VisGraphicsXform object. """
        self.translate_x = 0.0
        self.translate_y = 0.0
        self.translate_z = 0.0
        self.rotate_x = 0.0
        self.rotate_y = 0.0
        self.rotate_z = 0.0
        self.scale = 1.0

    def reset(self):
        """ Reset the transformation. """
        self.translate_x = 0.0
        self.translate_y = 0.0
        self.translate_z = 0.0
        self.rotate_x = 0.0
        self.rotate_y = 0.0
        self.rotate_z = 0.0
        self.scale = 1.0

    def set(self, xform):
        """ Set the transformation using a VisGraphicsXform object. """
        self.translate_x = xform.translate_x 
        self.translate_y = xform.translate_y 
        self.translate_z = xform.translate_z 
        self.rotate_x = xform.rotate_x
        self.rotate_y = xform.rotate_y
        self.rotate_z = xform.rotate_z
        self.scale = xform.scale 

class VisGraphicsPick(object):
    """ This class stores the picking data for the graphics scene. 

        Attributes:
            active (bool): If true then we are perfominf a pick. 
            color (List[Float]): The list of 4 (RGBA) color values used to display the pick line. 
            line_width (Float): The width used to display the pick line. 
            point1 (List[Float]): The first point defining the pick line.
            point2 (List[Float]): The second point defining the pick line.
            show_pick_line (bool): If true then show the 3D line passing though the picked point.
            intersect_point (List[Float]): The picked point closest to the viewer.
            intersect_points (List[List[Float]]): The list of picked points.

        Geometry in a scene can be selected by picking a point in the graphics window. This 2D point in screen
        coordinates is then transformed into a 3D line that is then used to intersect the geometry in the scene.
    """
    def __init__(self):
        self.active = False
        self.show_pick_line = False
        self.point1 = []
        self.point2 = []
        self.line_width = 2.0
        self.color = [0.5,1.0,0.5,1.0]
        self.intersect_point = None
        self.intersect_points = []

    def get_last_point(self):
        """ Get the last picked point closest to the viewer. """
        if not self.intersect_points:
            return None
        min_dist = 0
        min_point = None
        for i,point in enumerate(self.intersect_points): 
            v = [point[j] - self.point1[j] for j in xrange(03)]
            dist = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
            if (i == 0) or (dist < min_dist):
                min_dist = dist 
                min_point = point
        #__for point in self.intersect_points
        return min_point

    def render(self):
        """ Show the picked line and intersection point. """
        if self.show_pick_line and (len(self.point1) != 0):
            glLineWidth(self.line_width)
            glColor4fv(self.color)
            glBegin(GL_LINES)
            glVertex3dv(self.point1)
            glVertex3dv(self.point2)
            glEnd()

        if self.intersect_point: 
            self.draw_cross(self.intersect_point) 

    def draw_cross(self,point,width=0.08,color=[0.8,0.8,0.5,1.0],line_width=2.0):
        """ Draw a 3D cross at the given point. """
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

class VisGraphics(object):
    """ This class is used to manage a graphics scene, mouse events/actions and picking. 

        Attributes:
            action (MouseActions): The current mouse action. 
            center (List[Float]): The center of rotation. 
            cmd (VisComman): The command object used to write graphics actions to a file. 
            current_x (int): The current cursor x position.
            current_y (int): The current cursor y position.
            extent (VisExtent): The extent of the grahics scene. 
            height (int): The height of the graphics window.
            initial_xform (VisGraphicsXform): The transformation object storing the initial graphics scene rotation, 
                translation and scaling.
            menu (VisMenu): The menu object that manages the popup menu.
            pick (VisGraphicsPick): The pick object storing pick information.
            render_geometry (Dict[VisGeometry]): The list of geometry to render.
            title (String): The title of the graphics window.
            width (int): The width of the graphics window.
            xform (VisGraphicsXform): The transformation object storing the graphics scene rotation, translation and scaling.
            x_start (int): The x screen position at the start of a mouse move event. 
            y_start (int): The y screen position at the start of a mouse move event. 
    """
    def __init__(self, title, cmd):
        """ Initialize a VisGraphics object.

            Arguments:
                title (String): The graphics window title.
                cmd (VisCommand): The command processing object. 
        """
        self.title = title 
        self.command = cmd 
        self.dna_structure = None
        self.extent = VisExtent()
        self.pick = VisGraphicsPick()
        self.xform = VisGraphicsXform()
        self.initial_xform = VisGraphicsXform()
        self.center = [0.0,0.0,0.0] 
        self.width = 600
        self.height = 600
        self.action = MouseActions.NONE
        self.current_x = 0
        self.current_y = 0
        self.x_start = 0.0
        self.y_start = 0.0
        self.spectrum_colors = None
        self.menu = None 
        self.render_geometry = {}
        self._logger = logging.getLogger(__name__)

    def start_interactive(self):
        """ Start interactive graphics. """
        self._logger.info("Start graphics.")
        # Print help.
        self.help()
        # Turn the flow of control over to GLUT.
        glutMainLoop()

    def set_extent(self, extent):
        """ Set the graphics scene extent. """
        xmin,xmax,ymin,ymax,zmin,zmax = extent.get_bounds()
        self._logger.info("Extent: xmin %f  ymin %f  zmin %f" % (xmin, ymin, zmin))
        self._logger.info("        xmax %f  ymax %f  zmax %f" % (xmax, ymax, zmax))
        self.extent.set(xmin,xmax,ymin,ymax,zmin,zmax)
        cx,cy,cz = extent.get_center()
        self.center[0] = cx
        self.center[1] = cy
        self.center[2] = cz

    def set_center(self, point):
        """ Set the center of rotation. """
        cx,cy,cz = self.extent.get_center()
        self.center[0] = point[0] 
        self.center[1] = point[1] 
        self.center[2] = point[2] 
        self.xform.translate_x = cx - point[0] 
        self.xform.translate_y = cy - point[1] 
        self.xform.translate_z = cz - point[2] 
        self.xform.scale = 10.0
        self.command.generate_graphics_cmd("center", point)
        glutPostRedisplay()

    def center_on_pick(self):
        """ Set the center of rotation from a picked point. """
        pt = self.pick.get_last_point()
        if pt: 
            self.set_center(pt)
    
    def help(self):
        """ Print out help for mouse buttons and keyboard keys. """
        print("\n")
        print(" =================== Mouse buttons ===================")
        print(" Left mouse button - Rotate about the the view x-y axes ")
        print(" Middle mouse button - Translate the view in x-y ")
        print(" Right mouse button - Display the popup menu")
        print(" Shift-Left mouse button - Rotate about the view z axis ")
        print("\n")
        print(" ======================= Keys ========================")
        print(" esc - Quit ")
        print(" c - Query geometry and set center of rotation to selected point on the geometry")
        print(" q - Query geometry")
        print(" r - Reset view ")
        print(" + - Zoom in")
        print(" - - Zoom out")
        print("\n")
        VisMenu.help()

    def initialize_graphics(self):
        """ Initialize graphics. """
        # GLUT Window Initialization.
        glutInit()
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB| GLUT_DEPTH)

        # Set up window size and position.
        glutInitWindowSize (self.width, self.height)
        glutInitWindowPosition (0,0)
        #glutCreateWindow("Nanodesign Visualizer          File \'%s\' " % self.title)
        glutCreateWindow(self.title)

        # Initialize graphics view. 
        self.init_view()

        # Register callbacks.
        glutReshapeFunc(self.reshape)
        glutDisplayFunc(self.render)
        glutMouseFunc(self.mouse)
        glutMotionFunc(self.motion)
        glutPassiveMotionFunc(self.passive_motion)
        glutKeyboardFunc(self.keyboard)
        glutSpecialFunc(self.special_function)

    def add_render_geometry(self, geometry):
        """ Add a geometry to the render list. """
        self.render_geometry[geometry.id] = geometry

    def passive_motion(self, x, y):
        """ Process a passive mouse motion event. """
        self.current_x = x
        self.current_y = y

    def motion(self, x, y):
        """ Process a mouse motion event with a pressed button. """
        rs = 0.5
        dx = x - self.x_start
        dy = y - self.y_start

        if (self.action == MouseActions.ROTATE_XY):
            self.xform.rotate_x += rs*dy
            self.xform.rotate_y += rs*dx

        elif (self.action == MouseActions.ROTATE_Z):
            self.xform.rotate_z += rs*(dx+dy)

        elif (self.action == MouseActions.TRANSLATE_XY):
            wx, wy, wz = self.extent.get_widths()
            sx = wy / 200.0
            sy = wy / 200.0
            self.xform.translate_x +=  sx*dx
            self.xform.translate_y += -sy*dy

        elif (self.action == MouseActions.ZOOM):
            self.xform.scale += 0.05*(dx+dy)
            if ( self.xform.scale < 0.01):
                self.xform.scale = 0.01
        else:
            print("Unknown action\n", action)

        self.x_start = x
        self.y_start = y
        glutPostRedisplay()

    def special_function(self, key, x, y):
        """ Process a special function event: arrow keys."""
        if (key == GLUT_KEY_LEFT):
            self.xform.rotate_y += -90.0
        elif (key == GLUT_KEY_RIGHT):
            self.xform.rotate_y += 90.0
        elif (key == GLUT_KEY_UP):
            self.xform.rotate_x += -90.0
        elif (key == GLUT_KEY_DOWN):
            self.xform.rotate_x += 90.0
        glutPostRedisplay()

    def keyboard(self, key, x, y):
        """ Process a keyboard event. """

        # Set the center of rotation.
        if (key == 'c'):
            self.perform_pick(self.current_x, self.current_y)
            self.center_on_pick()

        # Query geometry.
        elif (key == 'q'):
            self.perform_pick(self.current_x, self.current_y)

        # Reset the view.
        elif (key == 'r'):
            self.reset_view()

        # Zoom in. 
        elif ((key == '+') or (key == 'i')):
            self.xform.scale += 0.1
            glutPostRedisplay()

        # Zoom out. 
        elif ((key == '-') or (key == 'o')):
            self.xform.scale -= 0.1
            glutPostRedisplay()
 
        # Esc key: quit.
        elif key == '\x1b':
            sys.exit(0)

    def mouse(self, button, state, x, y):
        """ Process a mouse button event. """
        if (button == GLUT_LEFT_BUTTON):
            if (glutGetModifiers() == GLUT_ACTIVE_SHIFT):
                self.action = MouseActions.ROTATE_Z
            elif (glutGetModifiers() == GLUT_ACTIVE_CTRL):
                self.action = MouseActions.ZOOM
            else:
                self.action = MouseActions.ROTATE_XY
        elif (button == GLUT_MIDDLE_BUTTON):
            self.action = MouseActions.TRANSLATE_XY
        elif (button == GLUT_RIGHT_BUTTON):
            self.action = MouseActions.ZOOM
        self.x_start = x
        self.y_start = y

    def clear_selections(self):
        """ Clear selections for all geometry. """
        for geom in self.render_geometry.values():
            geom.selected = False
        #__for geom in self.render_geometry.values()

    def perform_pick(self, x, y):
        """ Perform a pick operation. """
        self.clear_selections()
        self._logger.debug("Perform pick at screen (%d,%d) " % (x, y))

        # Convert the picked screen (x,y) into two 3D coordinates 
        # defining a line passing through the screen point.
        self.unproject_point(x, y)

        # Determine the geometry in the selected area by re-rendering it.
        glSelectBuffer(1000)
        glRenderMode(GL_SELECT)
        glInitNames()
        glPushName(0)
        self.pick.active = True
        self.display()

        # Process the list of geometry IDs determined at the picked point. 
        hit_buffer = glRenderMode(GL_RENDER)
        self._logger.debug("Hit record buffer size %d" % len(hit_buffer)) 
        num_hits = len(hit_buffer)
        picked_geoms = []
        for i,hit_record in enumerate(hit_buffer):
            min_depth, max_depth, names = hit_record
            self._logger.debug("Hit record min_depth %g max_depth %g   names %s" % (min_depth, max_depth, names)) 
            for id in names:
                if id in self.render_geometry:
                    geom = self.render_geometry[id]
                    picked_geoms.append(geom)
                    if geom.name:
                        name = geom.name
                    else:
                        name = " "
                    self._logger.debug("Selected geom id %s name %s " % (id, name))
            #__for id in names
        #__for hit_record in buffer:
        self.pick.active = False

        # Determine the intersection of the 3D line defined by the 
        # pick with the selected geometry.
        self.pick.intersect_points = []
        intersect_geoms = []
        self._logger.debug(" ")
        for geom in picked_geoms:
            self._logger.debug("Intersect selected geom name %s " % (geom.name))
            if geom.intersect_line(self.pick.point1, self.pick.point2):
                ipt = geom.intersect_point
                if not ipt:
                    continue
                self.pick.intersect_points.append(ipt)
                intersect_geoms.append(geom)
        #__for geom in picked_geoms
  
        # Find the closest point to the viewer.
        min_dist = None
        min_geom = None
        min_ipt = None
        for i,ipt in enumerate(self.pick.intersect_points):
            v = [self.pick.point1[j] - ipt[j] for j in xrange(0,3)]
            dist = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
            if (min_dist == None) or (dist < min_dist): 
                min_dist = dist
                min_geom = intersect_geoms[i]
                min_ipt = ipt
        #__for i,ipt in enumerate(self.pick.intersect_points)
 
        # Set the closest geometry to be selected.
        if (min_dist != None):
            self._logger.info("Selected geometry \'%s\'  point (%g, %g, %g)" % (min_geom.name, min_ipt[0], min_ipt[1], 
                 min_ipt[2]))
            min_geom.selected = True 
            if min_geom.selected_callback != None:
                min_geom.selected_callback(min_geom, min_geom.selected_entity)
            self.pick.intersect_point = min_ipt
        #__if (min_dist != None)

        # Re-display the scene.
        self.display()

    def unproject_point(self, sx, sy):
        """ Unproject a screen point. 

            Arguments:
                sx (int): The screen x coordinate. 
                sy (int): The screen y coordinate. 

            The two 3D end points of a line passing through (sx,sy) are calculated.
        """
        glLoadIdentity()
        glPushMatrix()
        cx = self.center[0]
        cy = self.center[1]
        cz = self.center[2]
        tx = cx + self.xform.translate_x
        ty = cy + self.xform.translate_y
        tz = cz + self.xform.translate_z
        glTranslatef(tx, ty, tz)
        glRotatef(self.xform.rotate_x, 1.0, 0.0, 0.0)
        glRotatef(self.xform.rotate_y, 0.0, 1.0, 0.0)
        glRotatef(self.xform.rotate_z, 0.0, 0.0, 1.0)
        glScalef(self.xform.scale, self.xform.scale, self.xform.scale)
        glTranslatef(-cx, -cy, -cz);

        # Get the viewing matrices.
        viewport = glGetIntegerv(GL_VIEWPORT)
        model_matrix = glGetDoublev(GL_MODELVIEW_MATRIX)
        proj_matrix = glGetDoublev(GL_PROJECTION_MATRIX)

        # Unproject the screen points.
        x = float(sx)
        y = float(viewport[3] - float(sy))
        z = 0.0
        wx1,wy1,wz1 = gluUnProject(x, y, z, model_matrix, proj_matrix, viewport)
        z = 1.0
        wx2,wy2,wz2 = gluUnProject(x, y, z, model_matrix, proj_matrix, viewport)

        # Set the endpoints for the picking line.
        #self.pick.show_pick_line = True
        self.pick.point1 = [wx1,wy1,wz1]
        self.pick.point2 = [wx2,wy2,wz2]
        glPopMatrix()

    def reset_view(self):
        """ Reset the view. """
        self.xform.set(self.initial_xform)
        cx,cy,cz = self.extent.get_center()
        self.center[0] = cx
        self.center[1] = cy
        self.center[2] = cz
        glutPostRedisplay()

    def reshape(self, width, height):
        """ Process a window reshape event. """
        self.height = height
        self.width = width
        self.display()
        # For the first reshape update the menu with selections from commands.
        # This is needed to make sure that the menus are fully initialized
        # and can then be modified. A flag in the 'menu' object makes sure
        # the update is done only once.
        self.menu.update()

    def display(self):
        """ Display the geometry defined for the graphics scene. """
        self.set_viewport(self.width, self.height)
        self.render()

    def set_viewport(self, width, height):
        """ Set the graphics viewport. """
        glViewport(0, 0, width, height)

        # Set up viewing transformation.
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        xmin,xmax,ymin,ymax,zmin,zmax = self.extent.get_bounds()
        cx, cy, cz = self.extent.get_center()
        dx, dy, dz = self.extent.get_widths()

        if (dx > dy):
            max_dim = dx
        else:
            max_dim = dy
        if (dz > max_dim):
            max_dim = dz

        sx = 1.0
        sy = 1.0

        if (width <= height):
            sy = float(height) / float(width)
            dy = dx * sy
            ymin = cy - dy / 2.0;
            ymax = cy + dy / 2.0;
        else:
            sx = float(width) / float(height)
            dx = dy * sx
            xmin = cx - dx / 2.0;
            xmax = cx + dx / 2.0;

        oxmin = cx - sx*max_dim;
        oxmax = cx + sx*max_dim;
        oymin = cy - sy*max_dim;
        oymax = cy + sy*max_dim;
        ozmin = cz - 100.0*max_dim;
        ozmax = cz + 100.0*max_dim;

        # If picking is active then set the picking matrix to define a picking window.
        if self.pick.active: 
            viewport = glGetIntegerv(GL_VIEWPORT)
            gluPickMatrix(float(self.current_x), float(viewport[3]-self.current_y), 5.0, 5.0, viewport)

        glOrtho(oxmin, oxmax, oymin, oymax, ozmin, ozmax)
        glMatrixMode(GL_MODELVIEW)

    def render(self):
        """ Render the geometry defined for the graphics scene. """
        # Clear frame buffer and depth buffer
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        #glClearColor(0.40, 0.60, 0.80, 1.0)
        glClearColor(0.75, 0.75, 0.75, 1.0)

        # Set transformations.
        glLoadIdentity()
        glPushMatrix ();
        cx = self.center[0]
        cy = self.center[1]
        cz = self.center[2]
        tx = cx + self.xform.translate_x
        ty = cy + self.xform.translate_y
        tz = cz + self.xform.translate_z
        glTranslatef(tx, ty, tz)
        glRotatef(self.xform.rotate_x, 1.0, 0.0, 0.0)
        glRotatef(self.xform.rotate_y, 0.0, 1.0, 0.0)
        glRotatef(self.xform.rotate_z, 0.0, 0.0, 1.0)
        glScalef(self.xform.scale, self.xform.scale, self.xform.scale)
        glTranslatef(-cx, -cy, -cz);

        # Render opaque geometry.
        pick_id = 1
        for geom in self.render_geometry.values():
            if not geom.transparent:
                if self.pick.active:
                    glLoadName(geom.id)
                geom.render()
                pick_id += 1
        #__for geom in self.render_geometry.values()

        # Render transparent geometry.
        for geom in self.render_geometry.values():
            if geom.transparent:
                if self.pick.active:
                    glLoadName(geom.id)
                geom.render()
                pick_id += 1
        #__for geom in self.render_geometry.values()

        glFlush()

        # Render picked point.
        self.pick.render()
        glPopMatrix();

        if not self.pick.active: 
            glutSwapBuffers()

    def init_view(self):
        """ Initialize lighting and rendering parameters. """
        ambient = [0.3, 0.3, 0.3, 0.0]
        glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
        glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, ambient);
        glEnable(GL_COLOR_MATERIAL);

        light_ambient =  [0.3, 0.3, 0.3, 0.0]
        light_diffuse =  [0.5, 0.5, 0.5, 1.0]
        light_specular =  [1.0, 1.0, 1.0, 1.0]
        light_position =  [0.0, 0.0, 1.0, 0.0]

        glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse)
        glLightfv(GL_LIGHT0, GL_POSITION, light_position)
        glEnable(GL_LIGHT0)

        glEnable (GL_NORMALIZE);
        glEnable (GL_LIGHTING);
        glShadeModel(GL_SMOOTH);

        glDepthFunc (GL_LESS);
        glEnable (GL_DEPTH_TEST);

        glEnable(GL_CULL_FACE);

        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glMatrixMode (GL_MODELVIEW);

        self.xform.set(self.initial_xform)
    #__def init_view

    def map_value_to_color(self, colors, vmin, vmax, value):
        num_colors = len(colors)
        dv = vmax - vmin
        if dv != 0.0: 
            f = (num_colors-1) / dv
            ci = int(f*(value - vmin))
            if ci < 0: ci = 0
            if ci > num_colors-1: ci = num_colors-1
            color = colors[ci]
        else:
            color = colors[0]
        #__if dv != 0.0

        return color[:]
    #__def get_color_map

    def get_spectrum_colors(self):
        """ Get a spectrum color map. """
        if self.spectrum_colors:
            return self.spectrum_colors

        num_colors = 16
        inc = 240.0 / (num_colors - 1)
        hue = 240.0
        self.spectrum_colors = []

        # Generate a list of RGB colors by incrementing hue.
        for i in xrange(0,num_colors):
            rgb = self.hsv_to_rgb(hue, 1.0, 1.0)
            self.spectrum_colors.append(rgb)
            hue -= inc
        #__for i in xrange(0,num_colors)

        return self.spectrum_colors
    #__def get_spectrum_colors 

    def hsv_to_rgb(self, h, s, v):
        """ Convert a hsv color to rgb. """
        if s == 0.0:
            return [v, v, v]
        else:
          if h == 360.0:
              h = 0.0;
          h = h / 60.0
          i = int(h)
          f = h - i
          p = v * (1.0 - s)
          q = v * (1.0 - (s*f))
          t = v * (1.0 - (s * (1.0 - f)))

          if i == 0:
              rgb = [v, t, p]
          elif i == 1:
              rgb = [q, v, p]
          elif i == 2:
              rgb = [p, v, t]
          elif i == 3:
              rgb = [p, q, v]
          elif i == 4:
              rgb = [t, p, v]
          elif i == 5:
              rgb = [v, p, q]
        #__if s == 0.0

        return rgb
    #__def hsv_to_rgb
      
#__class VisGraphics
