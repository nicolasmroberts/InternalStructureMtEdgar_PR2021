


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

# In statistics, 'inference' is a process of extrapolating from a data set to 
# the larger population that it represents. In particular, the mean of the 
# data is our best estimate of the population mean, but there is some 
# uncertainty in that estimate. Confidence regions and hypothesis tests are 
# two closely related ways of studying that uncertainty. This tutorial 
# develops confidence regions and hypothesis tests for lines.



### PROBLEM ###

# Load a data set of 23 foliation-lineation pairs from the western Idaho shear 
# zone (Giorgis and Tikoff, 2004). In this tutorial we'll treat only the 
# lineation directions. (A more detailed treatment can be found in the 
# orientation inference tutorial.)
wiszData <- geoDataFromFile("data/wiszFollins.tsv")
lineEqualAreaPlot(wiszData$direction)

# Giorgis and Tikoff (2004) explain the deformation in the western Idaho shear 
# zone as a shortening-dominated, homogeneous monoclinic transpression along a 
# NS-striking, vertical shear plane. Such a deformation predicts vertical 
# lineations.
wiszNull <- geoCartesianFromTrendPlungeDeg(c(0, 90))
lineEqualAreaPlotTwo(wiszData$direction, list(wiszNull), colorB="red")

# Of course, the observed lineations observed are not exactly vertical, just 
# because of random noise in the data. But is that the whole story? Can the 
# data be explained as random variation about an essentially vertical answer? 
# In other words, could wiszNull be the mean of the population, from which the 
# data arose?

# It is worth noting that all of our inference methods assume that the data 
# are 'independent and identically distributed' (IID). That is, they all arise 
# from the same underlying population, and they constitute a 'random' sample 
# from that population (as opposed to a systematically biased sample, for 
# example).



### WATSON DISTRIBUTION ###

# One approach is to assume that the data arose from a Watson distribution and 
# use Watson-specific methods. (This data set doesn't look very Watson-y, but 
# let's forge ahead anyway, because it's an instructive example.)
wiszWatson <- lineWatsonInference(wiszData$direction)
wiszWatson

# The following code produces a 'p-value', which is the probability of seeing 
# data like the western Idaho shear zone data (or data 'more extreme' than 
# them) if the mean lineation was vertical. If p < 0.05 (say) then we reject 
# the null hypothesis that the mean lineation is vertical.
wiszWatson$pvalue(wiszNull)

# Confidence regions are closely related to hypothesis tests. The 95% 
# confidence region consists of all null hypotheses that would NOT be rejected 
# at the p = 0.05 significance threshold. Here is that idea turned into code.
wiszConf <- lineWatson(mu=lineProjectedMean(wiszData$direction), kappa=10, n=10000)
wiszConf <- Filter(function(u) {wiszWatson$pvalue(u) >= 0.05}, wiszConf)
lineEqualAreaPlotTwo(wiszConf, list(wiszNull), colorB="red", shapeA=".")



### BINGHAM DISTRIBUTION ###

# The problem with the Watson treatment above is that the data might not have 
# arisen from a Watson distribution. In fact, the data cloud is elongated 
# enough that I doubt it. So let's try the more flexible Bingham distribution.
wiszBing <- lineBinghamInference(wiszData$direction, numPoints=100)
wiszMean <- lineProjectedMean(wiszData$direction)
lineEqualAreaPlotTwo(list(wiszMean), list(wiszNull), colorB="red", curves=list(wiszBing$points))

# The Bingham mean confidence region is reported as two angles, pointing from 
# the mean toward the other two principal directions from lineMeanScatter.
wiszBing$angles / degree

# We know that p < 0.05, because the null hypothesis is outside the 95% 
# confidence region. We don't know the p-value more precisely than that, 
# because this method does not report precise p-values.

# Unfortunately, there is little reason to believe that these data are Bingham-
# distributed. And n = 23 is a pretty small sample size for these methods. 
# (Last I checked, the Bingham functions in Allmendinger's and Cardozo's 
# Stereonet software balked for n < 25. And that might be wise.) So it would 
# be nice to have another approach.



### BOOTSTRAPPING ###

# In this section we try another approach called 'bootstrapping'. This method 
# is attractive because it is easy to implement and makes no distributional 
# assumptions.

# We resample the data set with replacement, to form a new data set of the 
# same size as the old one, but with some lines repeated and others omitted. 
# Intuitively, this new data set is like the original data set mutated a bit, 
# and its mean is like the original mean mutated a bit. We repeat this process 
# many times, to get a bunch of mutated means.
wiszBoot <- lineBootstrapInference(wiszData$direction, numBoots=10000, numPoints=100)
lineEqualAreaPlotTwo(wiszBoot$us, list(wiszNull), colorB="red")

# Bootstrapping theory (e.g., Efron and Tibshirani, 1993) says that these 
# mutated means quantify the uncertainty in the population mean (in a certain 
# way). We take the middle 95% of the bootstrapped means to form a 95% 
# confidence region for the population mean. Because the null hypothesis is 
# outside the 95% confidence region, we reject it as the population mean with 
# p < 0.05.
lineEqualAreaPlotTwo(list(lineProjectedMean(wiszData$direction)), list(wiszNull), colorB="red", curves=list(wiszBoot$points))

# To get the p-value more precisely, use this part of the output.
wiszBoot$pvalue(wiszNull)

# As an aside, let's use these results to illustrate some concepts. The 
# boundary of the 95% confidence region consists of exactly those null 
# hypotheses with p-value p == 0.05. The following code demonstrates this idea 
# by computing the p-value at a bunch of boundary points.
sapply(wiszBoot$points, wiszBoot$pvalue)

# Similarly, the points outside the region correspond to hypotheses that are 
# rejected (p < 0.05), and points inside the region correspond to hypotheses 
# that are not rejected (p > 0.05). The following code demonstrates this idea 
# by computing the p-value for each bootstrapped mean, and then computing 
# quantiles for those p-values. For example, 5% of the bootstrapped means 
# produce p < 0.05, and hence 5% of the bootstrapped means lie outside the 95% 
# confidence region, as we intended.
quantile(sapply(wiszBoot$us, wiszBoot$pvalue), probs=c(0.00, 0.05, 0.25, 0.50, 0.75, 1.00))



### WHAT DOES IT MEAN, GEOLOGICALLY? ###

# We now have multiple lines of argumentation suggesting that the western 
# Idaho shear zone lineations are not vertical. Therefore they are 
# incompatible with the proposed model. That model should be modified or 
# discarded.

# What specifically is wrong with the model? In other words, how should we 
# alter the model, so that it passes this statistical testing? Statistics does 
# not tell us that. That's the geologist's job.



### HOW TO REPORT AN INFERENCE ###

# When publishing inference results in a paper, do not merely report that you 
# reject or fail-to-reject the null hypothesis. Also report:
# 
# A. the specific method that you are using
# B. any important details --- for example, when bootstrapping, the number of 
#    bootstrap samples
# C. your chosen confidence level (usually 95%)
# D. the p-value (for methods that supply it)
# E. some idea of the effect size, such as a picture showing how far away from 
#    the confidence region the null hypothesis is.
# 
# The goal is to report enough information that the reader can determine 
# whether she agrees with your argument and reproduce your calculations if 
# necessary.



### OTHER OPTIONS ###

# The only problem with bootstrapping is that (like pretty much all 
# statistical methods) it doesn't work well on small data sets. The technique 
# of parametric bootstrapping can sometimes work around this problem.

# Another approach is Bayesian Markov chain Monte Carlo, which produces a 
# credible region (as opposed to a confidence region). I use it heavily in my 
# work, but we won't use it in these tutorials.



### CONCLUSION ###

# To draw inferences from your data about the larger population that it 
# represents is a little tricky. Some methods assume particular distributions 
# or large sample sizes. But you can use inference to corroborate or refute 
# geological explanations for data sets.


