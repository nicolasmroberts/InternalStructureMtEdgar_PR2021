


# Copyright 2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not 
# use this file except in compliance with the License. You may obtain a copy 
# of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software 
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
# License for the specific language governing permissions and limitations 
# under the License.



### INTRODUCTION ###

# This tutorial is an introduction to the basics of ellipsoids: How to 
# visualize them, how to describe their orientation, size, and shape, etc. We 
# introduce the crucial idea of degrees of freedom.



### IF YOU'VE JUST (RE)STARTED R ###

# If you've just restarted R, then you may need to set your working directory. 
# Your current working directory is named next to the word 'Console' at the 
# top of RStudio's Console pane. It should be the 'geologyGeometry' directory, 
# which contains subdirectories 'data', 'library', etc. If it isn't, then go 
# to the Files pane, navigate to the geologyGeometry directory, click the 
# 'More' menu, and choose 'Set As Working Directory'.

# Then click anywhere on the line of code below, and press RStudio's Run button 
# to execute it. What does it do? It loads the geologyGeometry library into R's 
# memory.
source("library/all.R")



### FINITE STRAIN ###

# Just for an example, let's make a finite strain ellipsoid by deforming a 
# spherical ball by a homogeneous monoclinic transpression (Fossen and Tikoff, 
# 1993).
fsePGT <- defMonoclinicPGT(gamma=2, logK=-1)
fseFinger <- fsePGT %*% t(fsePGT)
fseRotA <- ellRotationAFromTensor(solve(fseFinger))
ellEllipsoidPlot(rots=list(fseRotA$rotation), as=list(fseRotA$a))

# The 'a' field denotes the semi-axis lengths of the ellipsoid. Because this 
# finite strain ellipsoid has the same volume as the unit ball, the product of 
# the semi-axis lengths is 1.
fseRotA$a
prod(fseRotA$a)

# Geologists often prefer to use other ways of describing ellipsoid magnitude. 
# Most of them use the (natural) logarithms of the semi-axis lengths. They sum 
# to 0 when the ellipsoid has the same volume as the unit sphere.
fseLogA <- log(fseRotA$a)
fseLogA
sum(fseLogA)

# Here are some of the ways in which geologists describe magnitude. We'll 
# focus on (A) volume, which is always positive, (B) octahedral shear strain, 
# which is 0 for spheres and positive for other ellipsoids, and (C) Lode's 
# parameter, which is -1 for prolate spheroids, 1 for oblate spheroids, and 
# between -1 and 1 otherwise.
ellVolume(fseLogA)
ellOctahedralShearStrain(fseLogA)
ellLodeNu(fseLogA)
ellJelinekP(fseLogA)
ellFlinnK(fseLogA)

# Almost any three measures of magnitude will do. There are just a few 
# redundancies to avoid. For example, octahedral shear strain and Jelinek's P 
# obey a simple mathematical relationship and hence are redundant.
exp(sqrt(2) * ellOctahedralShearStrain(fseLogA))
ellJelinekP(fseLogA)

# In any event, you need three numbers to completely describe the magnitude of 
# an ellipsoid. We say that ellipsoid magnitudes have three 'degrees of 
# freedom'. Frequently, but not always, we ignore volume. Or rather we 
# normalize our ellipsoids to have the same volume as the unit ball. Then only 
# two degrees of freedom remain for magnitude.

# The Hsu-Nadai plot is a wedge-shaped plot for visualizing magnitude without 
# volume. Octahedral shear strain is the radial coordinate, increasing from 0 
# at the vertex. Lode's nu is approximately the angular coordinate, increasing 
# from -1 at the left edge to 1 at the right edge. (Warning: Actually nu is 
# not the angular coordinate. The plot is constructed in a different way.)
ellHsuNadaiPlot(list(fseLogA))

# Our R library also offers ellFlinnPlot, ellLogFlinnPlot, ellJelinekPlot, etc.

# The axes of an ellipsoid also have direction. This equal-area plot shows 
# short axis as a circle, the intermediate axis as a triangle, and the long 
# axis as a square.
ellEqualAreaPlot(rots=list(fseRotA$rotation), as=list(fseRotA$a))

# If you've absorbed the lessons of our orientation tutorials, then you might 
# agree that we should be viewing ellipsoid orientations holistically, as 
# orientations. They are subject to the 4-fold line-in-plane symmetry.
ellEqualVolumePlot(rots=list(fseRotA$rotation), as=list(fseRotA$a))

# Ellipsoid orientations, like all orientations, have three degrees of 
# freedom. Hence un-normalized ellipsoids have six degrees of freedom in 
# total, and volume-normalized ellipsoids have five degrees of freedom in 
# total.



### SHAPE PREFERRED ORIENTATION ###

# Shape preferred orientation (SPO) is an ellipsoid that summarizes a bunch of 
# ellipsoids in a rock. Despite its name, SPO is an ellipsoid, not an 
# orientation. Frequently it is constructed from elliptical sections of clasts 
# or grains on multiple outcrop faces. In this example, however, we obtain SPO 
# by averaging spinel clasts measured by X-ray computed tomography. The data 
# come from 31 field stations in New Caledonia (Titus et al., 2011; Chatzaras 
# et al., in prep). The file newcalOPXSpinelSPO.R loads these data into a list 
# ncSpinels.
source("data/newcalOPXSpinelSPO.R")

# Here are the magnitudes. The ellipsoids tend to be more prolate than oblate. 
ellHsuNadaiPlot(lapply(ncSpinels, function(stn) stn$spinel$logA))

# The orientations are all over the place.
ellEqualAreaPlot(rots=lapply(ncSpinels, function(stn) stn$spinel$rotation), 
                 as=lapply(ncSpinels, function(stn) stn$spinel$a))
ellEqualVolumePlot(rots=lapply(ncSpinels, function(stn) stn$spinel$rotation), 
                   as=lapply(ncSpinels, function(stn) stn$spinel$a))



### ANISOTROPY OF MAGNETIC SUSCEPTIBILITY ###

# When we apply a magnetic field to a rock that contains magnetic minerals, 
# those minerals affect the field. Anisotropy of magnetic susceptibility (AMS) 
# is an ellipsoid that quantifies this effect. In other words, AMS packages 
# information about the magnetic fabric of the rock. Here are 27 AMS 
# ellipsoids from Cyprus (Titus et al., in prep).
cyprusAMS <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv", doNormalize=FALSE)

# The ellipsoids are so close to spherical that we must zoom the Hsu-Nadai 
# plot. The second command here says 'show only the part of the plot where 
# octahedral shear strain is less than 0.03'.
ellHsuNadaiPlot(cyprusAMS$logA)
ellHsuNadaiPlot(cyprusAMS$logA, es=0.03)

# The orientations are a little better behaved than in the SPO example above, 
# but just a little.
ellEqualAreaPlot(cyprusAMS$rotation, cyprusAMS$a)
ellEqualVolumePlot(cyprusAMS$rotation, cyprusAMS$a)



### OTHER EXAMPLES ###

# When studying stress, it is sometimes useful to imagine the tractions 
# arising on all possible planes (that is, planes with all possible unit pole 
# vectors). When the eigenvalues of the stress tensor are all positive or all 
# negative, these tractions form a traction ellipsoid.

# Ellipsoid-like quantities also show up in permeability of fractures, thermal 
# conductivity, electrical conductivity, etc.



### CONCLUSION ###

# Ellipsoids have five (if normalized) or six (if not) degrees of freedom, 
# including three for orientation and two or three for magnitude. There are 
# various ways to measure and visualize these degrees of freedom.


