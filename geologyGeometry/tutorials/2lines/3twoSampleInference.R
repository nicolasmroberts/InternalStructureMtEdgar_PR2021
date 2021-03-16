


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

# Recall that inference is about extrapolating from a data set to the 
# population that it represents. Sometimes we want to perform an inference on 
# two data sets, to ascertain whether their populations could be identical. 
# This tutorial illustrates three approaches for lines: the test of Wellner 
# (1979), bootstrapping the mean, and two methods based on the Watson 
# distribution.



### PRELIMINARY WORK ###

# Load our dike data from two data files (Titus, pers. comm.). Glue the data 
# sets together. Just for fun, plot the data with Kamb contours (3-sigma, 
# 6-sigma, 9-sigma, 12-sigma, by default).
site229 <- geoDataFromFile("data/cyprusDikesSite229.tsv")
site230 <- geoDataFromFile("data/cyprusDikesSite230.tsv")
ourDikes <- rbind(site229, site230)
lineKambPlot(ourDikes$pole)

# Load some more dike data from an earlier map (Bear and Morel, 1960).
mapDikes <- geoDataFromFile("data/cyprusDikes_ic_three_blocks.tsv")
lineKambPlot(mapDikes$pole)

# We're wondering whether our data agree with the published data. So let's 
# inspect both data sets together. Notice that the map dikes (blue) seem a 
# little steeper than our dikes (red), but the two data sets overlap each 
# other pretty thoroughly.
lineEqualAreaPlotTwo(ourDikes$pole, mapDikes$pole, colorA="red", colorB="blue")

# Compute their means and inspect those. Yes, the map dike mean is a little 
# steeper than our dike mean. Of course, there is no reason to expect the two 
# means to be exactly identical.
ourMean <- lineProjectedMean(ourDikes$pole)
mapMean <- lineProjectedMean(mapDikes$pole)
lineEqualAreaPlotTwo(list(ourMean), list(mapMean), colorA="red", colorB="blue")

# If the two data sets agree, then we can combine them to make a larger data 
# set. Conversely, if the two data sets are incompatible somehow, then we 
# should not combine them. To address this issue, we try two approaches, both 
# based on simulation: the permutation test of Wellner (1979) and 
# bootstrapping.



### WELLNER (1979) TEST ###

# Wellner (1979) defined a statistic, denoted T, that measures how different 
# two sets of lines are. Let's compute it for these two data sets. The number 
# that comes out shouldn't give you any insight right now. But if you used 
# this statistic a lot, on many data sets, then you might start to get some 
# intuition for what it means. And its real purpose is to be used in the 
# permutation test below.
lineWellner(ourDikes$pole, mapDikes$pole)

# We have 31 dikes of our own and 38 dikes taken from a map. We are trying to 
# determine whether the distinction between 'our dikes' and 'map dikes' is a 
# meaningful concept. So we randomly reassign the dikes to these two groups. 
# That is, of the 69 total dikes, we arbitrarily pick 31 to call 'ours', and 
# we call the other 38 'from the map'. Then we recompute T. Intuitively, if 
# the new value of T is greater than the original value of T, then the new 
# distinction between 'ours' and 'from the map' is more meaningful than the 
# original distinction. But the new distinction is meaningless, so the 
# original distinction was probably meaningless too. If the new value of T is 
# less than the original value, then we have a little evidence that the 
# original distinction was meaningful.

# So in practice we compute a large number (say, 10,000) permutations of the 
# data, and count how many of them produce T greater than the original T. If 
# fewer than 5% of them produce greater T, then the original distinction 
# between 'ours' and 'from the map' is statistically significant at the 95% 
# confidence level. 
lineWellnerInference(ourDikes$pole, mapDikes$pole, 10000)

# When I did Wellner's test, I got p = 0.0016 based on 10,000 permutations and 
# p = 0.00193 based on 100,000 permutations. Assuming that your results are 
# similar (they'll be a little different every time this test is run), we have 
# strong evidence that our dikes and the map dikes come from different 
# populations. We'll discuss what this result means later.



### BOOTSTRAPPING ###

# Let's explore the same question using another approach. For each of the two 
# data sets, we bootstrap the mean, to get an idea of its uncertainty. (See 
# the one-sample inference tutorial for more explanation of bootstrapping.) If 
# the two clouds of bootstrapped means didn't overlap at all, then we could be 
# very confident that the means were different. But in this example the two 
# clouds do overlap a bit.
ourInf <- lineBootstrapInference(ourDikes$pole, numBoots=10000, numPoints=50)
mapInf <- lineBootstrapInference(mapDikes$pole, numBoots=10000, numPoints=50)
lineEqualAreaPlotTwo(ourInf$us, mapInf$us, colorA="red", colorB="blue")

# So we need to look at those bootstrapped means more closely. For each data 
# set, we construct a 95% confidence ellipse that contains the middle 95% of 
# those bootstrapped means. Because the two ellipses do not overlap, we can 
# conclude that the two populations are different.
lineEqualAreaPlotTwo(
  list(ourInf$center), list(mapInf$center), colorA="red", colorB="blue",
  curves=list(ourInf$points, mapInf$points))

# To illustrate p-values again, let's compute the greatest p-value attained by 
# any pole in the map confidence region, according to the notion of p-value 
# produced by our dikes. And vice-versa. Intuitively, each data set rejects 
# the other.
max(sapply(mapInf$points, ourInf$pvalue))
max(sapply(ourInf$points, mapInf$pvalue))



### WATSON ###

# Another approach is to assume that the data come from a particular 
# distribution --- in this case, the Watson --- and use methods specific to 
# that distribution. The geologyGeometry library offers two methods for multi-
# sample Watson inference.

# Warning: These data sets do not appear isotropic about their means. So they 
# might not be well-fit by the Watson distribution. So maybe we shouldn't 
# apply Watson techniques at all. But we'll go ahead anyway, just as an 
# illustration.

# One inference method assumes large sample size and tends to produce overly 
# conservative results, even when n = 3000. And our sample sizes are nowhere 
# near that big.
length(ourDikes$pole)
length(mapDikes$pole)

# But let's try it anyway. We get a p-value around 0.03, rejecting the null 
# hypothesis that the two distributions are equal. And, because the method is 
# overly conservative, we can trust that rejection.
lineLargeMultiSampleWatsonInference(list(ourDikes$pole, mapDikes$pole))

# The other inference method assumes tight concentration and tends to produce 
# good results as long as |kappa| >= 10. So let's check our concentration.
lineWatsonMLE(ourDikes$pole)
lineWatsonMLE(mapDikes$pole)

# Because the kappas are big enough, this method should work well. And it 
# gives a p-value around 0.001.
lineConcentratedMultiSampleWatsonInference(list(ourDikes$pole, mapDikes$pole))



### CONCLUSION ###

# The two dike populations look pretty similar, but multiple statistical tests 
# were able to detect a significant difference. Sometimes statistics is 
# sharper than your naked eye.

# What does it mean, geologically, that the dikes from the two studies come 
# from different populations? We were hoping to combine them into one larger 
# data set. Is that still a good idea? I'd say no. Data from the two data sets 
# are not comparable. Perhaps something was different in the methodology of 
# the two studies. For example, maybe they sampled from different locations, 
# and that difference turns out to matter. Or maybe the measurement apparatus 
# (pocket transit compass with human operator) was different in the two 
# studies.


