#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 10:07:07 2018

@author: snk

img2cup: convert an image into a 3D-printable cup

"""

import numpy as np
from math import *
from scipy import misc
from scipy import ndimage
from stl import mesh

# Parameters, change as needed
IMAGE = "escher.jpg"
MAX_AMP = 0.2      # Since final radius is 1, this value should be between 0 and 1
SOBEL = False      # Apply edge detection filter. In most cases not a good idea.
INVERT = False
N_MIRROR = 0       # Mirror N-times along the X axis to get a better aspect ratio. 

# Read and process image
img = misc.imread(IMAGE, flatten=True)
img /= img.max()
if SOBEL:
    sx = ndimage.sobel(img, axis=0, mode='constant')
    sy = ndimage.sobel(img, axis=1, mode='constant')
    img = np.hypot(sx, sy)

if INVERT:
    img = 1 - img
    
for i in range(N_MIRROR):
    img_flip = np.fliplr(img)
    img = np.append(img, img_flip, axis=1)    


shape = img.shape
print(img.shape)
verts = np.zeros([len(img.flatten()) + 2, 3], dtype="float64")
res = (2 * pi) / shape[1]


# Verts...
print("Vertices...")
zmax = res * (shape[0])# / shape[1]
verts[1] = np.array([0., 0, res * (shape[0] - 1)])
z = 0.
vidx = 2
for y in range(shape[0]):
    for x in range(shape[1]):
        z = (y / shape[0]) * zmax
        x_v = sin(res * x)
        y_v = cos(res * x)
        amp =  1 + img[y,x] * MAX_AMP
        verts[vidx] = np.array([x_v * amp, y_v * amp, z])
        vidx += 1
# Faces...
faces = []
print("Faces...")

# Upper Fan
up_start = shape[0] * shape[1] - shape[1] + 2
for i in range(shape[1]):
    faces.append([1, up_start + i + 1 if i != shape[1]-1 else up_start, up_start + i])
# Lower Fan
for i in range(shape[1]):
    faces.append([0, 2 + i, 2 + i + 1 if i != shape[1]-1 else 2])

# Cylinder connection
for y in range(shape[0] -1):
    for x in range(shape[1]):
        xp = x + 1 if x != shape[1]-1 else 0
        a = 2 + shape[1] * y + x
        b = 2 + shape[1] * y + xp
        c = 2 + shape[1] * (y + 1) + x
        d = 2 + shape[1] * (y + 1) + xp
        faces.append([a, c, d])
        faces.append([a, d, b])

# Create the mesh
print("Write STL...")

faces = np.array(faces)
cup = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        cup.vectors[i][j] = verts[f[j],:]

# Write the mesh to file "cube.stl"
cup.save('cup.stl')

