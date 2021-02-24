#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 13:39:16 2021

@author: phseal
"""

import hs_alkane.alkane as alk
import time
import numpy as np
import math as m


#########################
#  Testing timer module #
#########################
alk.timer_init() # initialise timer
time.sleep(1)

alk.cvar.timer_qtime = 10    # 10 seconds of allowed execution time
alk.cvar.timer_closetime = 1 #  5 seconds needed to stop safely

safe = alk.timer_check_continuation()
if (safe == 1):
    print("Safe to continue")


##########################
#  Testing random module #
##########################
seed = 234289
print("Seeding Fortran RNG with seed = ",seed)
alk.random_set_random_seed(seed)
xi = alk.random_uniform_random()
print("Random number xi = ",xi)


##############################
#  Testing quaternion module #
##############################
unit_vector1 = np.array([1.0,0.0,0.0],dtype=np.float64)
unit_vector2 = np.array([0.0,0.0,1.0],dtype=np.float64)

# Find quaternion which rotates unit_vector1 onto unit_vector2
quat1 = alk.quat_get_minimum_arc(unit_vector1, unit_vector2)
print(quat1)

# Find quaternion which rotates by -pi/2 around y axis
axis = np.array([0.0,1.0,0.0],dtype=np.float64)
angle = -0.5*m.pi
quat2 = alk.quat_axis_angle_to_quat(axis, angle)
print(quat2)

# Find product of quat1 and quat2
quat3 = alk.quat_product(quat1, quat2, 1)
print(quat3)

# Find inverse of quat2
quat4 = alk.quat_inverse(quat2)

# Rotate unit_vector1 using quat1
unit_vector3 = alk.quat_conjugate_q_with_v(quat1, unit_vector1)
print(unit_vector3)

# Rotate unit_vector1 using quat2
unit_vector3 = alk.quat_conjugate_q_with_v(quat2, unit_vector1)
print(unit_vector3)

# Rotate back again using quat4
unit_vector1 = alk.quat_conjugate_q_with_v(quat4, unit_vector3)
print(unit_vector1)


#######################
#  Testing box module #
#######################
num_replicas = alk.box_get_num_boxes()
print("Number of replicas : ",num_replicas)

#alk.box_set_num_boxes(2)
#num_replicas = alk.box_get_num_boxes()[0]
#print("Number of replicas : ",num_replicas)
#alk.box_set_num_boxes(1)
#num_replicas = alk.box_get_num_boxes()[0]
#print("Number of replicas : ",num_replicas)

alk.box_initialise()

print(alk.box_get_cell(1))

cell_matrix = np.array([[ 7.120546, 0.000000,  0.000000],
                        [ 1.584668, 7.859930,  0.000000],
                        [-1.815922, 2.549568,-10.042803]], dtype=np.float64)


alk.box_set_cell(1, cell_matrix)
alk.box_update_recipmatrix(1)


#cell_matrix = np.array([[11,12,13],[21,22,23],[31,32,33]], dtype=np.float64)
alk.box_set_cell(1, cell_matrix)
print(alk.box_get_cell(1))

pos1 = np.array([0.2,0.2,0.2],dtype=np.float64)
pos2 = np.array([0.8,0.8,0.8],dtype=np.float64)

pos1c = alk.box_frac_to_cart(1, pos1)
pos2c = alk.box_frac_to_cart(1, pos2)

print("pos1 converted to Cartesian coords : ",pos1c)
print("pos2 converted to Cartesian coords : ",pos2c)

pos1 = alk.box_cart_to_frac(1,pos1c)
pos2 = alk.box_cart_to_frac(1,pos2c)

print("pos1c converted to Fractional coords : ",pos1)
print("pos2c converted to Fractional coords : ",pos2)

minv = alk.box_minimum_image(1, pos1c, pos2c)
print("Minimum image vector between pos1 and pos2 :",minv)

alk.box_set_pbc(0)
minv = alk.box_minimum_image(1, pos1c, pos2c)
print("Direct vector between pos1 and pos2        :",minv)
alk.box_set_pbc(1)

vol = alk.box_compute_volume(1)
print("Volume of simulation box : ",vol)

alk.box_set_bypass_link_cells(0)
alk.box_construct_link_cells(ibox=1)

alk.box_set_pbc(0)

############################
#  Testing alkane module  #
###########################
nchains = alk.alkane_get_nchains()
print("Number of chains : ",nchains)

nchains = 1
alk.alkane_set_nchains(nchains)
nchains = alk.alkane_get_nchains()
print("Number of chains : ",nchains)

nbeads = alk.alkane_get_nbeads()
print("Number of beads per chain : ",nbeads)
nbeads = 4
alk.alkane_set_nbeads(nbeads)
nbeads = alk.alkane_get_nbeads()
print("Number of beads per chain : ",nbeads)

alk.alkane_init()  # Inconsistent init vs initialise in box

# Grow a chain from scratch using CBMC
ifail = 1
while ifail != 0:
    rbfactor, ifail = alk.alkane_grow_chain(1, 1, new_conf=1)

print("Rosenbluth factor : ",rbfactor)
print("ifail : ",ifail)

# Build initial linked lists for this initial config
alk.alkane_construct_linked_lists(ibox=1)

# Numpy array which refers to chain coordinates
mychain = alk.alkane_get_chain(1,1)
print(mychain)
print()

orig_chain = mychain.copy()

# Trial move - translation
backup_chain = mychain.copy()
Pacc = alk.alkane_translate_chain(1,1)
print("Acceptance probability for move : ", Pacc)
print(mychain)
print()
print(mychain-backup_chain)

# Trial move - rotation
backup_chain = mychain.copy()
Pacc, quaternion = alk.alkane_rotate_chain(1,1,bond=1)
print("Acceptance probability for move : ", Pacc)
print(mychain)
print()
print(mychain-backup_chain)

# Trial move - torsion
backup_chain = mychain.copy()
Pacc, idihedral, angle = alk.alkane_bond_rotate(1,1,allow_flip=1)
print("Rotated dihedral bond no. ",idihedral," by ", angle, "radians")
print("Acceptance probability for move : ", Pacc)
print(mychain)
print()

# Check for chain overlaps
overlap = alk.alkane_check_chain_overlap(1)
if overlap!=0:
    print("Overlap between chains!")


# Check chain geometry
violated = alk.alkane_check_chain_geometry(1,1)
if violated != 0:
    print ("Geometry of chain ",1," in box ",1," violates constraints")

mychain[1][1] += 0.1

violated = alk.alkane_check_chain_geometry(1,1)
if violated != 0:
    print ("Geometry of chain ",1," in box ",1," violates constraints")

mychain[1][1] -= 0.1

violated = alk.alkane_check_chain_geometry(1,1)
if violated != 0:
    print ("Geometry of chain ",1," in box ",1," violates constraints")

# Test link cell update
for ibead, bead in enumerate(mychain):
    print("Updating linked lists for bead ",ibead)
    alk.alkane_update_linked_lists(ibead+1,1,1,orig_chain[ibead],bead)
    
alk.alkane_destroy()
