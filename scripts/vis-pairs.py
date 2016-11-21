#!/usr/bin/env python

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

""" This script is used to interactively visualize a SimDNA pairs file using OpenGL.

    Strands are displayed in different colors with C5' atoms displayed as spheres. Base pairings are
    displayed as lines beteween C5' atoms.

    New geometry can be rendered by adding functions to the VisGraphics.render() method.

    The script uses the PyOpenGL package for a Python binding to OpenGL.
"""
import os
import sys
import random 
import numpy as np
from numpy import linalg

try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    from OpenGL.GLUT import *

except ImportError as e:
    print "Could not import PyOpenGL."
    raise e

class Base(object):
    """ This class stores information for a base in a strand.  """
    def __init__(self, strand_id, base_id, coord, paired_base_id, paired_strand_id):
        self.strand_id = strand_id
        self.base_id = base_id
        self.coord = [coord[0], coord[1], coord[2]]
        self.paired_base_id = paired_base_id
        self.paired_strand_id = paired_strand_id

class DnaModel(object):
    """ This class stores information for a SimDNA model read from a .paris file. """
    def __init__(self):
        self.num_bases = 0
        self.pairs = []
        self.strands = {}
        self.strand_colors = {}

    def read_pairs_file(self,file_name):
        """ Read in a SimDNA pairs file. """
        print(">>> read pairs file %s" %  file_name) 
        with open(file_name, 'r') as infile:
            lines = infile.readlines()
        self.num_bases = int(lines[0])
        print(">>> number of bases %d" %  self.num_bases) 

        for i in xrange(1,len(lines)):
            line = lines[i].strip('\n')
            line = ' '.join(line.split())
            tokens = line.split(" ")
            n = 0
            strand_id = int(tokens[n])
            n += 1
            base_id = int(tokens[n])
            n += 1
            coord = [float(tokens[n]), float(tokens[n+1]), float(tokens[n+2])]
            n += 3
            paired_strand_id = int(tokens[n])
            n += 1
            paired_base_id = int(tokens[n])
            base = Base(strand_id, base_id, coord, paired_base_id, paired_strand_id)
            self.pairs.append(base)
        #__for i in xrange(1,len(lines))

    def get_strands(self):
        """ Get a list of strands. """
        if not self.strands:
            random.seed(1.0)
            for base in self.pairs:
                strand_id = base.strand_id
                if strand_id not in self.strands:
                    self.strands[strand_id] = []
                    color = (random.random(), random.random(), random.random())
                    self.strand_colors[strand_id] = color
                self.strands[strand_id].append(base)
            print(">>> number of strands %d" % len(self.strands))
        #__if not self.strands
        return  self.strands

#__class DnaModel

#====================================================================================#
#                               g r a p h i c s                                      #
#====================================================================================#
class MouseActions:
   NONE         = -1
   ROTATE_XY    = 0
   ROTATE_Z     = 1
   TRANSLATE_XY = 2
   TRANSLATE_Z  = 3
   ZOOM         = 4

class VisGraphicsXform(object):
    """ This class stores the transformation for the graphics scene. """
    def __init__(self):
        self.translate_x = 0.0
        self.translate_y = 0.0
        self.translate_z = 0.0
        self.rotate_x = 0.0
        self.rotate_y = 0.0
        self.rotate_z = 0.0
        self.scale = 1.0

    def reset(self):
        self.translate_x = 0.0
        self.translate_y = 0.0
        self.translate_z = 0.0
        self.rotate_x = 0.0
        self.rotate_y = 0.0
        self.rotate_z = 0.0
        self.scale = 1.0

class VisGraphicsExtent(object):
    """ This class stores extent information. """
    def __init__(self):
        self.is_set = False
        self.xmin = 0.0
        self.ymin = 0.0
        self.zmin = 0.0
        self.xmax = 1.0
        self.ymax = 1.0
        self.zmax = 1.0
        self.cx = 0.0
        self.cy = 0.0
        self.cz = 0.0

    def set(self,x,y,z):
        self.xmin = x
        self.ymin = y
        self.zmin = z
        self.xmax = x
        self.ymax = y
        self.zmax = z
        self.is_set = True

    def reset(self):
        self.xmin = 1e6
        self.ymin = 1e6
        self.zmin = 1e6
        self.xmax = -1e6
        self.ymax = -1e6
        self.zmax = -1e6

    def update(self,x,y,z):
        if not self.is_set:
            self.set(x,y,z)
            return

        if (x < self.xmin):
            self.xmin = x
        elif (x > self.xmax):
            self.xmax = x

        if (y < self.ymin):
            self.ymin = y
        elif (y > self.ymax):
            self.ymax = y

        if (z < self.zmin):
            self.zmin = z
        elif (z > self.zmax):
            self.zmax = z

    def get_bounds(self):
        return self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax

    def get_center(self):
        return (self.xmin+self.xmax)/2.0, (self.ymin+self.ymax)/2.0, (self.zmin+self.zmax)/2.0

    def get_widths(self):
        return (self.xmax-self.xmin), (self.ymax-self.ymin), (self.zmax-self.zmin)

class VisGraphicsPairs(object):
    def __init__(self, colors1, pair1_coords, colors2, pair2_coords, mid_coords):
       self.colors1 = colors1
       self.pair1_coords = pair1_coords
       self.colors2 = colors2 
       self.pair2_coords = pair2_coords 
       self.mid_coords = mid_coords 

class VisGraphics(object):
    """ This class stores information for a graphics scene and renders a SimDNA model. """ 
    def __init__(self, model):
        self.dna_model = model
        self.extent = VisGraphicsExtent()
        self.xform = VisGraphicsXform()
        self.width = 600
        self.height = 600
        self.action = MouseActions.NONE
        self.x_start = 0.0 
        self.y_start = 0.0
        self.pairs_geometry = None

    def start_interactive(self):
        print(">>> start graphics") 
        self.set_extent()
        self.initialize_graphics()

    def set_extent(self):
        for base in self.dna_model.pairs: 
            self.extent.update(base.coord[0], base.coord[1], base.coord[2])
        xmin,xmax,ymin,ymax,zmin,zmax = self.extent.get_bounds()
        print("[VisGraphics] extent: xmin %f  ymin %f  zmin %f" % (xmin, ymin, zmin))
        print("                      xmax %f  ymax %f  zmax %f" % (xmax, ymax, zmax))

    def initialize_graphics(self):
        # GLUT Window Initialization.
        glutInit()
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB| GLUT_DEPTH)

        # Set up window size and position.
        glutInitWindowSize (self.width, self.height)
        glutInitWindowPosition (0,0)
        glutCreateWindow ("SimDNA pairs visualizer")

        # Initialize graphics view. 
        self.init_view()

        # Register callbacks.
        glutReshapeFunc(self.reshape)
        glutDisplayFunc(self.render)
        glutMouseFunc(self.mouse)
        glutMotionFunc(self.motion)
        glutKeyboardFunc(self.keyboard)

        # Turn the flow of control over to GLUT.
        glutMainLoop()

    def motion(self, x, y):
        """ Process a mouse motion event. """
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
            print("unknown action\n", action)

        self.x_start = x
        self.y_start = y
        glutPostRedisplay()

    def keyboard(self, key, x, y):
        """ Process a keyboard event. """

        if (key == 'r'):
            self.reset_view()

        if ((key == '+') or (key == 'i')):
            self.xform.scale += 0.1
            glutPostRedisplay()

        if ((key == '-') or (key == 'o')):
            self.xform.scale -= 0.1
            glutPostRedisplay()

        if (key == 'q'):
            exit(0)

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

    def reset_view(self):
        """ Reset the view. """
        self.xform.reset()
        glutPostRedisplay()

    def reshape(self, width, height):
        """ Process a window reshape event. """
        self.height = height
        self.width = width
        #print("[reshape] width %d  height %d " % (width, height))
        self.display()

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

        glOrtho(oxmin, oxmax, oymin, oymax, ozmin, ozmax)
        glMatrixMode(GL_MODELVIEW)

    def render(self):
        """ Render the geometry defined for the graphics scene. """
        # Clear frame buffer and depth buffer
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glClearColor (0.7, 0.8, 0.9, 1.0)

        glLoadIdentity()
        glPushMatrix ();
        cx, cy, cz = self.extent.get_center()
        tx = cx + self.xform.translate_x
        ty = cy + self.xform.translate_y
        tz = cz + self.xform.translate_z
        glTranslatef(tx, ty, tz)

        glRotatef(self.xform.rotate_x, 1.0, 0.0, 0.0)
        glRotatef(self.xform.rotate_y, 0.0, 1.0, 0.0)
        glRotatef(self.xform.rotate_z, 0.0, 0.0, 1.0)
    
        glScalef(self.xform.scale, self.xform.scale, self.xform.scale)

        glTranslatef(-cx, -cy, -cz);

        # Render geometry.
        self.render_strands()
        self.render_pairs()

        glPopMatrix ();
        glutSwapBuffers()

    def render_strands(self):
        """ Render the SimDNA model strands. """
        strands = self.dna_model.get_strands()

        for id,base_list in strands.items(): 
            glLineWidth(4.0)
            glBegin(GL_LINE_STRIP)
            glColor3fv(self.dna_model.strand_colors[id])
            for base in base_list: 
                glVertex3f(base.coord[0], base.coord[1], base.coord[2])
            glEnd()
            for base in base_list:
                self.render_sphere(base.coord, self.dna_model.strand_colors[id])

    def render_pairs(self):
        """ Render the SimDNA model paired bases. """

        # Create pairs geometry as two lines between pairs colored
        # by the colors of the strands they join.

        if not self.pairs_geometry: 
            self.pairs_geometry = [] 
            strands = self.dna_model.get_strands()
            pair_coords1 = {}
            pair_coords2 = {}

            for strand_id,base_list in strands.items(): 
                for base in base_list: 
                    paired_base_id = base.paired_base_id 
                    paired_strand_id = base.paired_strand_id 

                    if (paired_base_id != -1):
                        paired_strand = strands[paired_strand_id]
                        paired_base = paired_strand[paired_base_id-1] 
                        x = (base.coord[0] + paired_base.coord[0]) / 2.0
                        y = (base.coord[1] + paired_base.coord[1]) / 2.0 
                        z = (base.coord[2] + paired_base.coord[2]) / 2.0

                        color1 = self.dna_model.strand_colors[strand_id]
                        if color1 not in pair_coords1: 
                            pair_coords1[color1] = []
                        pair_coords1[color1].append( (base.coord[0], base.coord[1], base.coord[2], x, y, z) )

                        color2 = self.dna_model.strand_colors[paired_strand_id]
                        if color2 not in pair_coords2: 
                            pair_coords2[color2] = []
                        pair_coords2[color2].append( (paired_base.coord[0], paired_base.coord[1], paired_base.coord[2], x, y, z) )
                    #__if (paired_base_id != -1)
                #__for base in base_list
            #__for strand_id,base_list in strands.items()

            self.pairs_geometry = [pair_coords1, pair_coords2]
        #__if not self.pairs_geometry: 

        glLineWidth(2.0)
        glBegin(GL_LINES)
        pair_coords1, pair_coords2 = self.pairs_geometry
        for color, pair_coords in pair_coords1.items():
            glColor3fv(color)
            for pair_coord in pair_coords:
                coord = pair_coord[0:3]
                mid_coord = pair_coord[3:]
                glVertex3fv(coord)
                glVertex3fv(mid_coord)
        glEnd()

    def render_sphere(self, center, color):
        """ Render a sphere. """
        radius = 1.0
        sides = 8
        glPushMatrix()
        try:
            glColor3fv(color)
            glTranslatef(center[0],center[1],center[2])
            glutSolidSphere(radius, sides, sides)

        finally:
            glPopMatrix()

    def init_view(self):
        """ Initialize lighting and rendering parameters. """ 
        ambient = [0.3, 0.3, 0.3, 0.0]
        glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
        glLightModeli (GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
        glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, ambient);
        glEnable (GL_COLOR_MATERIAL);

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

def help():
    print(" =================== mouse buttons ===================")
    print(" Left mouse button - rotate about model x-y axes ")
    print(" Middle mouse button - translate in x-y ")
    print(" Right mouse button - zoom ")
    print(" Shift-Left mouse button - rotate about model z axis ")
    print(" =================== keys ===================")
    print(" q - quit ")
    print(" r - reset view ")
    print(" + - zoom in")
    print(" - - zoom out")

def main():

    if len(sys.argv) != 2:
        print("**** ERROR: wrong number of arguments")
        print("\nUsage: vis-pairs.py <filename>")
        sys.exit(0)

    file_name = sys.argv[1]
    dna_model = DnaModel()
    dna_model.read_pairs_file(file_name)

    help()

    vis_gr = VisGraphics(dna_model)
    vis_gr.start_interactive()

if __name__ == "__main__":
    main()
    
