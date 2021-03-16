


# Copyright 2016 Joshua R. Davis
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

# Many geological data types are not simply lines (or rays), but pairs (or 
# triples) of lines. One example is foliation-lineation pairs. The foliation 
# can be described by its pole line, the lineation can also be described as a 
# line, and those two lines are inherently perpendicular. Another example is 
# axial planes and hinge lines of folds. As this tutorial demonstrates, we 
# should NOT treat the lines within each pair separately, because that would 
# ignore the perpendicularity.



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



### MEAN PLANE POLE VS. MEAN LINE ###

# We're going to randomly generate a synthetic data set of 10 data, where each 
# datum is a plane containing a line (like a foliation-lineation). Then we'll 
# compute the mean plane pole and the mean line direction and the angle 
# between them. It should be 90 degrees...?
rotations <- rotIsotropicFisher(rotUniform(), 3, 10)
poles <- lapply(rotations, function(r) r[1,])
directions <- lapply(rotations, function(r) r[2,])
lineEqualAreaPlotTwo(poles, directions, shapeA="c", shapeB="s")
acos(dot(lineProjectedMean(poles), lineProjectedMean(directions))) / degree

# With your mouse, select all five of the above lines of code at once. Copy 
# them, paste them into the Console pane, and press Return to execute them. 
# Then press the up-arrow key and Return, to execute the entire code block 
# again. Each time you will get a slightly different angle.

# This code block runs 1,000 such experiments, making a histogram of the 
# angles that result. Typically there are many results more than 5 degrees 
# away from perpendicular. Sometimes there are even some results 15 degrees 
# away.
perpendicularMeanTrial <- function(concen, n) {
  rotations <- rotIsotropicFisher(rotUniform(), concen, n)
  poles <- lapply(rotations, function(r) r[1,])
  directions <- lapply(rotations, function(r) r[2,])
  acos(dot(lineProjectedMean(poles), lineProjectedMean(directions))) / degree
}
hist(replicate(1000, perpendicularMeanTrial(3, 10)))

# Computing the sample mean is a really basic operation, that might affect a 
# variety of more sophisticated techniques (such as bootstrapping). So we want 
# our notion of sample mean to be well-behaved, with reliable mathematical 
# properties. Computing the mean directions of planes and lines separately 
# isn't going to be acceptable.



### HIDDEN OUTLIER EXAMPLE ###

# Here's another synthetic data set of plane-line pairs --- 28 of them. The 
# plane poles are plotted as circles and the lines are plotted as squares. I 
# claim that there is a dramatic outlier in this data set --- one plane-line 
# pair that is quite different from the others. Can you see it in this plot? 
# (No, you can't.)
outlierData <- geoDataFromFile("data/synthFoldsOutlier.csv")
lineEqualAreaPlotTwo(outlierData$pole, outlierData$direction, shapeA="c", shapeB="s")

# The problem with that plot is that it divorces each plane from its 
# associated line. So a lot of information is lost. Here is another plot, that 
# shows each line on its plane. In theory, you should be able to see the 
# outlier now. In practice, it's still difficult to see.
lineEqualAreaPlot(outlierData$direction, curves=lapply(outlierData$pole, rayGreatCircle, 72), shapes="s")

# Now we arbitrarily color some of the plane-line pairs blue, some of them 
# green, and one of them red. Do you see why the red datum is an outlier? 
# (We'll see it even more clearly in the next few sections.)
lineEqualAreaPlot(c(outlierData$pole, outlierData$direction), colors=as.character(outlierData$color),
                  shapes=c(replicate(nrow(outlierData), "c"), replicate(nrow(outlierData), "s")))



### CONCLUSION ###

# If your data consist of multiple inter-dependent lines or rays, then you 
# should not treat them as if they were independent.


