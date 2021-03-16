


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

# In the preceding tutorial we learned the 'right' way to compute the mean and 
# dispersion of an orientational data set. But tools for orientations are so 
# scarce that we must sometimes resort to a less-than-ideal way.



### THE WRONG WRONG WAY ###

# Load a set foliation-lineation orientations from the western Idaho shear 
# zone (Giorgis and Tikoff, 2004). Also compute and illustrate the Frechet 
# mean. Because of the 4-fold symmetry, the everything appears four times in a 
# rotation plot.
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
wiszMeanVar <- oriMeanVariance(wiszData$rotation, group=oriLineInPlaneGroup)
wiszCurves <- lapply(
  oriNearestRepresentatives(wiszData$rotation, wiszMeanVar$mean, group=oriLineInPlaneGroup),
  rotGeodesicPoints, wiszMeanVar$mean)
oriEqualAnglePlot(wiszData$rotation, curves=wiszCurves, group=oriLineInPlaneGroup)

# What if we weren't being careful, and we treated our orientations as mere 
# rotations? First, here's the plot. The data set looks bimodal, even though 
# it's unimodal.
rotEqualAnglePlot(wiszData$rotation)

# Here's the rotational mean as minimizing the size of a spider in that plot.
wiszMeanVarRot <- rotMeanVariance(wiszData$rotation)
wiszCurvesRot <- lapply(wiszData$rotation, rotGeodesicPoints, wiszMeanVarRot$mean)
rotEqualAnglePlot(wiszData$rotation, curves=wiszCurvesRot)

# The spider looks much bigger than in the orientational treatment. The 
# variances confirm it.
wiszMeanVar$variance
wiszMeanVarRot$variance



### THE RIGHT WRONG WAY ###

# Remembering now that we should be dealing with orientations, let's view all 
# four symmetric copies of that spider. This plot reveals the problem: By 
# ignoring symmetry, the rotation-only treatment 'crosses between symmetric 
# copies' of the data when it shouldn't.
oriEqualAnglePlot(wiszData$rotation, curves=wiszCurvesRot, group=oriLineInPlaneGroup)

# In this example, we can take another tack: At the start, ensure that all of 
# the orientations are represented by rotations in a single symmetric copy. 
# Then treat them as rotations.
wiszRots <- oriNearestRepresentatives(wiszData$rotation, wiszData$rotation[[1]], group=oriLineInPlaneGroup)
wiszMeanVarOkay <- rotMeanVariance(wiszRots)
wiszCurvesOkay <- lapply(
  wiszRots,
  rotGeodesicPoints, wiszMeanVarOkay$mean)
rotEqualAnglePlot(wiszRots, curves=wiszCurvesOkay)

# When we symmetrize that plot, we get the same plot as we did in the true, 
# orientational treatment.
oriEqualAnglePlot(wiszData$rotation, curves=wiszCurvesOkay, group=oriLineInPlaneGroup)

# And the variances confirm that we're getting the right answer.
wiszMeanVar$variance
wiszMeanVarOkay$variance

# The lesson here is: If your data are tightly concentrated enough, then 
# isolating the symmetric copies is not difficult, and working within one 
# symmetric copy often yields good results.

# This strategy is common in electron backscatter diffraction (EBSD), for 
# example. Let's take another look at the Moine thrust intra-grain quartz 
# orientations of (Strine and Wojtal, 2004; Michels et al., 2015). Six tight, 
# beautiful symmetric copies.
michelsData <- geoDataFromFile("data/moine_one_grainABCxyz.tsv")
oriEqualVolumePlot(michelsData$rotation, group=oriTrigonalTrapezohedralGroup, simplePoints=TRUE)

# But structural orientations are often so dispersed that we can't easily 
# discern the symmetric copies. Here are some faults with slip from Cyprus 
# (Davis and Titus, 2017). Where are the two symmetric copies?
slickData <- geoDataFromFile("data/cyprusSlicks2008005.tsv")
oriEqualVolumePlot(slickData$rotation, group=oriRayInPlaneGroup)

# So the other part of the lesson is: Structural geologists cannot always 
# ignore symmetry in orientational data. Our default approach should be one 
# that handles symmetry. Only sometimes can we cheat on the symmetry.



### SO LET'S CHEAT THEN ###

# Let's return to our western Idaho shear zone foliation-lineation data set.
oriEqualVolumePlot(wiszRots, group=oriLineInPlaneGroup)

# Remember that we've pre-processed the data to choose representative 
# rotations lying in one symmetric copy. So intrinsic methods will give the 
# same results for these rotations as for the corresponding orientations.
rotEqualVolumePlot(wiszRots)

# Maximum likelihood estimation (MLE) of the matrix Fisher distribution 
# parameters is a method for rotations. And it's not intrinsic, so we're not 
# allowed to use it on orientations. But these orientations are concentrated 
# enough that treating them as rotations is approximately okay. So let's do 
# the MLE. Among other things, it gives a concentration matrix (kHat) whose 
# eigenvalues quantify the dispersion in the data set.
wiszFisher <- rotFisherMLE(wiszRots)
eigen(wiszFisher$kHat, symmetric=TRUE)$values

# In the next tutorial we will learn why this specific calculation is useful.



### SOMETHING THAT RESEMBLES CHEATING ###

# When your data are tightly concentrated, you can approximate them as points 
# in the tangent space at the mean. Then principal component analysis (PCA) in 
# that tangent space gives you yet another measure of anisotropic dispersion.
wiszPCA <- rotLeftPrincipalComponentAnalysis(wiszRots, wiszMeanVarOkay$mean, numPoints=5)
wiszPCA$magnitudes
rotEqualAnglePlot(wiszRots, curves=wiszPCA$curves, simplePoints=TRUE)

# Theoretical aside: We're using rotation methods, but this PCA concept is 
# actually intrinsic to the geometry of the space of orientations, so it is a 
# legitimate orientation technique. We're not actually cheating.

# Anyway, here is the symmetrized version of that last plot.
oriEqualAnglePlot(wiszRots, group=oriLineInPlaneGroup, curves=wiszPCA$curves, simplePoints=TRUE)



### CONCLUSION ###

# For tightly concentrated orientation data, intrinsic orientation methods 
# work exactly as intrinsic rotation methods do. Sometimes we cheat and use 
# non-intrinsic rotation methods on orientations too.


