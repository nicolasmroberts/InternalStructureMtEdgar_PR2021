


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

# Now that we have a system for converting plane-line pairs into rotations (up 
# to 4-fold symmetry), we can plot those rotation using the equal-angle and 
# equal-volume rotation plots. They can be regarded as equal-angle and equal-
# area hemispherical plots, 'one dimension up'. Their only major drawback is 
# that each datum appears four times because of the symmetry.



### WESTERN IDAHO SHEAR ZONE EXAMPLE ###

# Load the 23 foliation-lineation data from the western Idaho shear zone 
# (Giorgis and Tikoff, 2004). For each datum there are four representative 
# rotations. We plot all 4 * 23 = 92 rotations in an equal-angle rotation plot 
# and an equal-volume rotation plot. You might want to view them side by side.
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
oriEqualAnglePlot(wiszData$rotation, group=oriLineInPlaneGroup)
oriEqualVolumePlot(wiszData$rotation, group=oriLineInPlaneGroup)

# For comparison, here are the foliation poles and lineation directions 
# separately on an equal-area hemispherical plot.
lineEqualAreaPlotTwo(wiszData$pole, wiszData$direction, shapeA="c", shapeB="s")

# If you're going to ask a question about directions (either the lineations or 
# the foliation poles), then you use the equal-area plot. These directions are 
# very difficult to read from the rotation plots.

# But if you're going to ask a question about orientations (the foliation-
# lineations holistically), then you use the rotation plots, because the 
# equal-area plot is a much less faithful view of the space of orientations.

# I'm not proposing that you stop using equal-area plots. I'm encouraging you 
# to pick up the rotation plots, as a complementary tool in your tool set.



### HIDDEN OUTLIER EXAMPLE ###

# Here's the hidden outlier example again.
outlierData <- geoDataFromFile("data/synthFoldsOutlier.csv")
lineEqualAreaPlotTwo(outlierData$pole, outlierData$direction, shapeA="c", shapeB="s")

# Each datum plots four times, producing four symmetric copies of the data. In 
# each copy, the outlier is unmistakable.
oriEqualAnglePlot(outlierData$rotation, group=oriLineInPlaneGroup)
oriEqualVolumePlot(outlierData$rotation, group=oriLineInPlaneGroup)



### NOT EASY? ###

# In the rest of these tutorials, when we need to plot rotational or 
# orientational data, we will usually employ the equal-angle or equal-volume 
# plot. For example, the equal-volume plot will help us determine whether a 
# data set is unimodal (one region of high density, centered on a mean) or 
# bimodal (two regions of high density, so that the mean doesn't make sense). 
# The equal-volume plot depicts distance and shape relationships imperfectly, 
# but much better than an equal-area plot of planes and lines separately, for 
# example.

# With all that said, you may still feel that equal-angle and equal-volume 
# plots of rotational versions of structural data is all too strange. Here is 
# my advice:

# A. Equal-angle and equal-area hemispherical plots were not obvious either, 
# when you first saw them. They required practice. So do these plots. And they 
# require more practice than the equal-area plot, because orientations are 
# inherently more complicated than rays or lines.

# B. Even if you were in love with the equal-angle and equal-volume plots, I 
# would not recommend that you use them exclusively. Always view your data in 
# multiple kinds of plot, to minimize your chances of overlooking important 
# features of your data set.

# C. Even if you've already vowed that you'll never use any of these plots 
# again, the rest of these tutorials still offer a variety of useful 
# statistical methods that have little to do with plotting.



### CONCLUSION ###

# Depending on the question you're asking, you might want to use the equal-
# area (or equal-angle) hemispherical plot, or the equal-angle or equal-volume 
# rotation plot. But if you're studying orientations holistically, pick one 
# of the latter two.


