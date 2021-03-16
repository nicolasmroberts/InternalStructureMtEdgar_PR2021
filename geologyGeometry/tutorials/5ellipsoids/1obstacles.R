


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

# In this tutorial we start thinking about analyzing ellipsoids. Various 
# obstacles arise. For example, how can we do statistics on ellipsoid 
# magnitudes? Simple calculations don't work well, fundamentally because 
# magnitudes don't form a vector space. How can we do statistics on ellipsoid 
# orientations? We have lots of tools for that, but those tools don't work 
# well when spheroids (or near-spheroids) are present in the data.

# We would really like to treat ellipsoids holistically, as unified objects, 
# rather than piecemeal. Fortunately, there is a way to treat ellipsoids 
# holistically as elements in a vector space. This idea forms the basis for 
# all of the methods in subsequent sections.



### COMPUTING WITH MAGNITUDES ###

# Here are the semi-axis lengths of two ellipsoids. Let's average them. R 
# gives an answer, but it's a bad answer. Why?
obstacleA <- c(5, 1, 2)
obstacleB <- c(1, 5, 2)
arithmeticMean(list(obstacleA, obstacleB))

# The easiest solution is to insist that the semi-axis lengths be recorded in 
# a certain order (ascending or descending).
sort(obstacleA)
sort(obstacleB)
arithmeticMean(list(sort(obstacleA), sort(obstacleB)))

# Geologists do this all the time. Frequently they go on to focus on the 
# longest semi-axis or the shortest semi-axis. So suppose that we have a bunch 
# of ellipsoids, and here are their longest semi-axis lengths:
obstacleA1s <- c(0.4, 0.5, 2.1, 0.3, 0.6, 0.1, 0.2)

# Using a common R function, here is the 95% confidence interval for the mean 
# of the longest axis. Why is this answer bad?
t.test(obstacleA1s)

# Statisticians run into this kind of problem a lot. If your data are positive 
# numbers, in a context where zero and negative numbers don't even make sense, 
# then you should a method that's aware of that fact. Frequently the solution 
# is to (A) transform to the logarithms of the numbers, (B) do statistics 
# there, and (C) transform the results back.
obstacleLogA1s <- log(obstacleA1s)
obstacleLogA1sTest <- t.test(obstacleLogA1s)
exp(obstacleLogA1sTest$estimate)
exp(as.numeric(obstacleLogA1sTest$conf.int))

# Mathematically, the fundamental reason why this tactic works is that the 
# logarithms live in a 'vector space'. A vector space is a setting where basic 
# arithmetic operations (addition, subtraction, scaling) make sense. Hence the 
# statistical calculations built from those operations also make sense.

# So to analyze entire size-shapes of ellipsoids, should we just take the logs 
# of their semi-axis lengths, and work with those vectors of logs? Sorry, but 
# no. These log-vectors are not just any 3D vectors. They have to be ordered 
# in descending (or ascending order). Consequently the set of all log-vectors 
# does not form a vector space. Statistics is not going to work well.



### COMPUTING WITH ORIENTATIONS ###

# Load a data set of 27 AMS ellipsoids from Cyprus (Titus et al., in prep). 
# Their orientations are quite widely dispersed.
cyprusAMS <- geoEllipsoidDataFromIRMFile("data/cyprus_AMS_groupF.tsv", doNormalize=FALSE)
ellEqualVolumePlot(cyprusAMS$rotation, cyprusAMS$a)
ellEqualAreaPlot(cyprusAMS$rotation, cyprusAMS$a)

# To some degree, the axes tend to align in three dense regions of the equal-
# area plot. You can see what I mean by forgetting which axis is which and 
# then Kamb contouring. Each dense region contains a mixture of short, 
# intermediate, and long axes.
lineKambPlot(unlist(lapply(cyprusAMS$rotation, function(r) list(r[1,], r[2,], r[3,])), recursive=FALSE))

# It's almost as if the short, intermediate, and long axes are getting mixed 
# up. This plot suggests a reason.
ellHsuNadaiPlot(cyprusAMS$logA)

# If an ellipsoid is clearly triaxial, meaning that its semi-axis lengths are 
# quite different, then its orientation is well-defined (up to 4-fold line-in-
# plane symmetry).

# But if an ellipsoid is spheroidal, meaning that two of its semi-axes have 
# the same length, then its orientation is ill-defined. A circle's worth of 3D 
# orientations, or equivalently a single 3D direction, describes the 
# spheroid's orientation.

# And when ellipsoids are close to spheroidal, then tiny errors in their 
# magnitudes can produce huge changes in their orientations. 'Mixing up' the 
# short, intermediate, and long axes is just what you expect to see in such a 
# data set.

# Therefore, when a data set contains a mixture of triaxial ellipsoids and 
# spheroids, finding a coherent way to analyze their orientations separately 
# from their magnitudes is difficult.



### THE EASIEST WAY TO OVERCOME THESE OBSTACLES ###

# In summary, it's not obvious how to do statistics on ellipsoid magnitudes. 
# And it's not obvious how to do statistics on ellipsoid orientations, when 
# their magnitudes are not cooperating.

# In a fantasy world, we would analyze ellipsoids holistically, in a way that 
# accounts for magnitude and orientation together. And that holistic treatment 
# would happen in a vector space, where statistics is easy.

# Remarkably, this world exists. There is a way of packaging an ellipsoid's 
# five or six degrees of freedom (depending on whether it is normalized or 
# not) into a 5D or 6D vector, which we call a 'log-ellipsoid vector'. Under 
# this system, each ellipsoid corresponds to one and only one vector, and each 
# vector corresponds to one and only one ellipsoid. The one-to-one 
# correspondence is quite well-behaved; for example, ellipsoids that are 
# similar correspond to vectors that are close together.

# This system facilitates a simple workflow for doing statistics with 
# ellipsoids:
# 
# A. Convert the ellipsoids into 5D or 6D vectors. In our library, you usually 
#    do this by asking the ellipsoid for its $vector.
# B. Do any kind of statistics on those vectors that you want. R has hundreds 
#    of things for you to try.
# C. Convert the vectors back into ellipsoids. In our library, this is the 
#    function ellEllipsoidFromVector, although sometimes more complicated 
#    operations are required.
# D. Inspect the ellipsoids to understand the result.
# 
# We'll see our first example in the next section.



### CONCLUSION ###

# There is a way to express ellipsoids as 5D or 6D vectors. That is the format 
# in which pretty much all ellipsoid statistics should be done.


