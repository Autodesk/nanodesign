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
This module is contains various functions used to create a DNA structure from a caDNAno DNA origami design file.
"""
import numpy as np
from math import sqrt,cos,acos,sin,asin,pi
from .common import CadnanoLatticeType

def generate_coordinates(dna_parameters, lattice_type, row, col, helix_num, scaffold_bases, staple_bases):
    """ Generate the axis coordinates, axis reference frames, and nucleotide coordinate for the
        virtual helix. 

        Arguments:
            dna_parameters (DnaParameters): The DNA parameters to use when creating the 3D geometry for the design.
            lattice_type (CadnanoLatticeType): The lattice type for this design.
            row (int): The caDNAno row number. 
            col (int): The caDNAno column number. 
            helix_num (int): The caDNAno virtual helix number. 
            scaffold_bases (List[DnaBase]): The list of scaffold bases for this helix. 
            staple_bases (List[DnaBase]): The list of staple bases for this helix. 

        Returns:
            axis_coords (NumPy Nx3 ndarray[float]): The coordinates of base nodes along the helix axis.
            axis_frames (NumPy 3x3xN ndarray[float]): The coordinate frames of base nodes alonge the helix axis.
            scaffold_coords ((NumPy Nx3 ndarray[float]): The coordinates of nucleotides along the DNA helix. 
            staple_coords ((NumPy Nx3 ndarray[float]): The coordinate of nucleotides along the DNA helix.

        The helix coordinates and reference frames are generated for all of the virtural helix
        positions that contain a base. The coordinates and refereance frames are also set for 
        scaffold and staple bases. 
    """
    r_strand = dna_parameters.helix_distance / 2.0 # half the distance between the axes of adjacent DNA helices.
    r_helix = dna_parameters.helix_radius          # radius of DNA helices (nm)
    dist_bp = dna_parameters.base_pair_rise        # rise between two neighboring base-pairs (nm)
    ang_bp = dna_parameters.base_pair_twist_angle  # twisting angle between two neighboring base-pairs (degrees)
    ang_minor = dna_parameters.minor_groove_angle  # angle of the minor groove (degrees)

    # Positions of the scaffold nucleotide and staple nucleotide
    # in the local reference frame.
    scaf_local = r_helix * np.array([cos(deg2rad(180-ang_minor/2)), sin(deg2rad(180-ang_minor/2)), 0.0]).transpose()
    stap_local = r_helix * np.array([cos(deg2rad(180+ang_minor/2)), sin(deg2rad(180+ang_minor/2)), 0.0]).transpose()

    # Set the helix start coordinates. 
    init_coord,init_ang = get_start_coordinates_angle(dna_parameters, lattice_type, row, col, helix_num)

    # Set the direction of the helix axis.
    if (helix_num % 2 == 0):
        e3 = np.array([0, 1, 0],dtype=float)
    else:
        e3 = np.array([0, -1, 0],dtype=float)

    # Create a list of sorted base positions.
    base_positions = set()
    for base in scaffold_bases: 
        base_positions.add(base.p)
    for base in staple_bases: 
        base_positions.add(base.p)

    # Compute helix axis coordinates and frames.
    num_base_positions = len(base_positions)
    axis_coords = np.zeros((num_base_positions,3), dtype=float)
    axis_frames = np.zeros((3,3,num_base_positions), dtype=float)
    pos_map = {}
    for i,p in enumerate(sorted(base_positions)): 
        pos_map[p] = i
        ref_coord = init_coord + np.array([0, dist_bp*p, 0])
        angle = init_ang + ang_bp*p
        # Base coordinate. 
        axis_coords[i,:] = ref_coord
        # Base orientation.
        e2 = np.array([cos(-deg2rad(angle)), 0, sin(-deg2rad(angle))])
        e1 = np.cross(e2, e3);
        axis_frames[:,:,i] = np.array([e1, e2, e3]).transpose();
    #__for base in staple_bases

    # Compute scaffold nucleotide positions and set base coordinates and frame.
    num_scaffold_bases = len(scaffold_bases)
    scaffold_coords = np.zeros((num_scaffold_bases,3), dtype=float)
    for i,base in enumerate(scaffold_bases): 
        j = pos_map[base.p]
        base.coordinates = axis_coords[j]
        base.ref_frame = axis_frames[:,:,j]
        scaffold_coords[i,:] = axis_coords[j,:] + np.dot(axis_frames[:,:,j], scaf_local)
        base.nt_coords = scaffold_coords[i]
    #__for base in scaffold_bases 

    # Compute staple nucleotide positions and set base coordinates and frame.
    num_staple_bases = len(staple_bases)
    staple_coords = np.zeros((num_staple_bases,3), dtype=float)
    for i,base in enumerate(staple_bases): 
        j = pos_map[base.p]
        base.coordinates = axis_coords[j]
        base.ref_frame = axis_frames[:,:,j]
        staple_coords[i,:] = axis_coords[j,:] + np.dot(axis_frames[:,:,j], stap_local)
        base.nt_coords = staple_coords[i]
    #__for base in stape_bases

    return axis_coords, axis_frames, scaffold_coords, staple_coords 
#__def generate_coordinates

def get_start_coordinates_angle(dna_parameters, lattice_type, row, col, helix_num):
    """ Get the start axis coordinates and angle for a virtual helix. 

        Arguments:
            dna_parameters (DnaParameters): The DNA parameters to use when creating the 3D geometry for the design.
            lattice_type (CadnanoLatticeType): The lattice type for this design.
            row (int): The caDNAno row number. 
            col (int): The caDNAno column number. 
            helix_num (int): The caDNAno virtual helix number. 
    """
    r_strand = dna_parameters.helix_distance / 2.0 # half the distance between the axes of adjacent DNA helices.
    dist_bp = dna_parameters.base_pair_rise        # rise between two neighboring base-pairs (nm)
    ang_bp = dna_parameters.base_pair_twist_angle  # twisting angle between two neighboring base-pairs (degrees)

    # Set the helix start coordinates, based on helix (row,col), and frame orientation, 
    # based on helix number (polarity).
    if lattice_type == CadnanoLatticeType.honeycomb:
        xpos = sqrt(3.0) * col * r_strand
        zpos = -3.0 * row * r_strand;
        if ( ((row % 2 == 0) and (col % 2 == 0)) or ((row % 2 != 0) and (col % 2 != 0))):
            zpos = zpos + r_strand
        if ( helix_num % 2 == 0 ):
            init_ang = -30 + ang_bp / 2
        else:
            init_ang = 150 + ang_bp / 2
    elif lattice_type == CadnanoLatticeType.square: 
        xpos =  2.0 * col * r_strand
        zpos = -2.0 * row * r_strand
        if helix_num % 2 == 0:
            init_ang = 180 + ang_bp / 2
        else:
            init_ang = 0 + ang_bp / 2
    #__if lattice_type == CadnanoLatticeType.honeycomb

    return np.array([xpos, 0.0, zpos]), init_ang
#__def get_start_coordinates_angle

def deg2rad(deg):
    """Convert degrees into radians. """
    rad = (pi/180)* deg;
    return rad
#__def deg2rad

def find_row(neigh, curr_bases):
    tmp = curr_bases - neigh
    s = np.sum(np.abs(tmp),1)
    ind = np.where(s==0)[0]
    return ind
#__def find_row

def vrrotmat2vec(R):
    """ Extract the equivalent rotation about an axis from a rotation matrix. """
    m00 = R[0,0]
    m01 = R[0,1]
    m02 = R[0,2]
    m10 = R[1,0]
    m11 = R[1,1]
    m12 = R[1,2]
    m20 = R[2,0]
    m21 = R[2,1]
    m22 = R[2,2]
    angle = acos(( m00 + m11 + m22 - 1)/2.0)
    x = (m21 - m12) / sqrt( pow(m21-m12,2) + pow(m02-m20,2) + pow(m10-m01,2) )
    y = (m02 - m20) / sqrt( pow(m21-m12,2) + pow(m02-m20,2) + pow(m10-m01,2) )
    z = (m10 - m01) / sqrt( pow(m21-m12,2) + pow(m02-m20,2) + pow(m10-m01,2) )
    return np.array([x,y,z],dtype=float),angle
#__def vrrotmat2vec

def vrrotvec2mat(axis, theta):
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
#__def vrrotvec2mat

def bp_interp(dnode_1, triad_1, dnode_2, triad_2, n):
    """ Interpolate the position (dnode) and orientation (triad) between two base pairs.

        Solve for rotation matrix R defined as the solution of: 
            R * triad_1 = triad_2 (multiply both sides by transpose(triad_1).
    """
    dnode_interp = np.zeros((n,3),dtype=float)
    triad_interp = np.zeros((3,3,n),dtype=float)

    # Solve for rotation matrix R.
    R = np.dot(triad_2,triad_1.T)

    # Get the equivalent rotation about an axis for R.
    a,theta = vrrotmat2vec(R)

    # Calculate for dNode_interp and triad_interp.
    for i in xrange(0,n):
        dnode_interp[i,:] = (dnode_1*(n+1-(i+1)) + dnode_2*(i+1)) / (n+1)
        angle = theta*(i+1)/(n+1)
        rot_mat = vrrotvec2mat(a, angle)
        triad_interp[:,:,i] = np.dot(rot_mat,triad_1)

    return dnode_interp, triad_interp

#__def bp_interp
