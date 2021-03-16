


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

# This tutorial begins our study of rotations in three-dimensional space. One 
# way to conceptualize a rotation is as an axis and an amount of rotation 
# about that axis. The axis is a ray, and positive rotation amounts occur in a 
# 'right-handed' sense about that ray by convention. If we further express the 
# axis ray using trend and plunge, then we have three angles --- trend, 
# plunge, amount --- describing the rotation.

# There are many other ways of recording and calculating with rotations. All 
# of them use at least three numbers. We say that rotations have three 
# 'degrees of freedom' or that the space of rotations is three-dimensional. 
# It follows that a faithful plot of rotations must be a 3D plot.

# Rotations are sometimes directly useful in geology. This tutorial gives a 
# couple of examples. But their main purpose is to form the theoretical 
# foundation for orientations.



### IF YOU'VE JUST (RE)STARTED R ###

# If you've just restarted R, then you may need to set your working directory. 
# Your current working directory is named next to the word 'Console' at the 
# top of RStudio's Console pane. It should be the 'geologyGeometry' directory, 
# which contains subdirectories 'data', 'library', etc. If it isn't, then go 
# to the Files pane, navigate to the geologyGeometry directory, click the 
# 'More' menu, and choose 'Set As Working Directory'.

# Then click anywhere on the line of code below, and press RStudio's Run button 
# to execute it. It loads the geologyGeometry library into R's memory.
source("library/all.R")



### PLATE MOTION RECONSTRUCTIONS ###

# The relative movement of two tectonic plates is rotational (rather than 
# translational). The axis of the rotation is called the Euler pole. 
# Engebretson et al. (1984) published reconstructions of relative plate 
# motions. We use their Farallon-Pacific relative motion as processed by 
# Prentice (1987).

# The farpacTs variable contains 18 dates, from 0 to 165 (in millions of years 
# before the present). The farpacRs variable contains 18 rotations, each 
# rendered as a 3x3 rotation matrix.
source("data/farpacRotations.R")
farpacTs
farpacRs

# I don't want to dwell on those rotation matrices, but let's just inspect the 
# first datum. It says that 0 million years ago the net rotation was about the 
# ray <0, 0, -1> through an angle of 0 radians. That's the trivial rotation 
# --- the rotation that doesn't actually rotate. From this example we can infer 
# that rotations are being measured relative to the present.
farpacTs[[1]]
rotAxisAngleFromMatrix(farpacRs[[1]])

# The fact that the rotations get larger as we go back in time suggests that 
# they are also being measured cumulatively (which is true).



### ANGLE-AXIS PLOT AND TWO VARIANTS ###

# Suppose that we're talking about a rotation with axis U (a ray or unit 
# vector) and amount A (an angle in radians). Consider the scaled vector A U. 
# Its direction is the same as U --- that is, the axis of rotation --- and its 
# length is A, the amount of rotation. We plot this vector A U as a point in 
# the usual way: by placing its tail at the origin and marking the point where 
# its head lands. In this way, every rotation plots as a point in the solid 
# ball of radius pi. This system for plotting rotations is called the 
# 'axis-angle plot'. Points on the boundary of the ball are antipodally 
# identified, just like the points on the boundary of a hemispherical plot.

# Here's the axis-angle plot of the Farallon-Pacific data. We're using color 
# to depict the time variable, from red (present) to magenta (distant past). 
# By playing with the plot window, you can learn how to rotate and zoom the 
# image.
rotAxisAnglePlot(farpacRs, colors=hues(farpacTs), boundaryAlpha=0.3)

# A question for you: Where is the pure-red dot and why?

# It turns out that the axis-angle plot does not enjoy an equal-angle or equal-
# volume property. However, it can be distorted into plots with those 
# properties. Try viewing the plots below side-by-side with the axis-angle 
# plot, to see how subtle the difference is. (It's like the subtle difference 
# between the equal-angle and equal-area hemispherical plots.)
rotEqualAnglePlot(farpacRs, colors=hues(farpacTs), boundaryAlpha=0.3)
rotEqualVolumePlot(farpacRs, colors=hues(farpacTs), boundaryAlpha=0.3)

# For comparison, here is a plot showing only the axis of rotation. It does 
# not show any information about the amount of rotation. Of course, it could 
# be annotated with the amounts, but even then it would be misleading. To see 
# what I mean, inspect two reddest dots.
rayEqualAreaPlot(lapply(farpacRs, function(r) rotAxisAngleFromMatrix(r)[1:3]), colors=hues(farpacTs))

# The two reddest dots are far apart in the equal-area plot, although the 
# equal-volume plot reveals that they should be close together. Viewing only 
# one aspect of the data --- in this case, the axis without the amount --- can 
# mislead us.



### INFERRING ROTATION FROM PALEOMAGNETIC DIRECTION ###

# Suppose that we've got a paleomagnetic direction with declination 308 
# degrees and inclination 21 degrees.
pmag <- geoCartesianFromTrendPlungeDeg(c(308, 21))
rayEqualAreaPlot(list(pmag))

# But the expected direction, for rocks of this age and latitude, is 
# declination 0 degrees and inclination 5 degrees, suggesting that the rocks 
# have been deformed.
should <- geoCartesianFromTrendPlungeDeg(c(0, 5))
rayEqualAreaPlot(list(should))

# Based on these two rays, can we infer the deformation? No. Deformations are 
# too rich a class to be constrained so easily. The answer would be highly non-
# unique.

# Here's a less ambitious question: If we assume that the deformation was 
# actually just a rigid rotation, then can we infer that rotation? No. The 
# problem is still slightly underconstrained. We must supply at least one 
# additional piece of information --- new data, an assumption, etc.

# Here's one idea: Of all the possible rotations, find the smallest one. By 
# reading and running this code, can you guess what each part does?
smallRot <- rotSmallestRotationFromTwoRays(should, pmag)
axisAngle <- rotAxisAngleFromMatrix(smallRot)
geoTrendPlungeDegFromCartesian(axisAngle[1:3])
axisAngle[[4]] / degree
rayEqualAreaPlot(list(should, pmag, axisAngle[1:3]), curves=list(rayGeodesicPoints(should, pmag)))
rayEqualAreaPlot(list(should, pmag, axisAngle[1:3]), curves=list(rayGreatCircle(axisAngle[1:3])))

# Another idea is: Suppose that we also have bedding data. The younged pole to 
# bedding is trend 218 degrees, plunge -14 degrees. But the expected younged 
# pole would have plunge -90 degrees (i.e. point straight up). The code here 
# gets a bit complicated, but don't worry about the details.
young <- geoCartesianFromTrendPlungeDeg(c(218, -14))
init <- rotProjectedMatrix(cbind(should, c(0, 0, 1), cross(should, c(0, 0, 1))))
final <- rotProjectedMatrix(cbind(pmag, young, cross(pmag, young)))
extraRot <- final %*% t(init)
axisAngle <- rotAxisAngleFromMatrix(extraRot)
geoTrendPlungeDegFromCartesian(axisAngle[1:3])
axisAngle[[4]] / degree
rayEqualAreaPlot(list(should, pmag, axisAngle[1:3], c(0, 0, 1), young), 
                 curves=list(raySmallCircle(axisAngle[1:3], rayDistance(axisAngle[1:3], pmag)), 
                             raySmallCircle(axisAngle[1:3], rayDistance(axisAngle[1:3], young))))

# Here are both of those deduced rotations in an equal-volume plot. Can you 
# tell which is which?
rotEqualVolumePlot(list(smallRot, extraRot))



### CONCLUSION ###

# Rotations have three degrees of freedom and thus benefit from 3D plots. In 
# particular, the equal-angle and equal-volume plots are much like the equal-
# angle and equal-area hemispherical plots, but in one higher dimension.


