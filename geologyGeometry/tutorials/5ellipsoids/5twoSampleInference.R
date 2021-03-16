


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

# Frequently we have two data sets concerning a single geologic system, and we 
# want to know whether those two data sets are telling us similar or disparate 
# information. In particular, it is common to ask whether the populations, 
# from which the data arise, have the same mean. To answer such questions we 
# need two-sample inference.



### EXPLORATION ###

# Load two X-ray computed tomography data sets from the Huemul plutonic 
# complex, Chile (Garibaldi et al., in prep). The two data sets measure mafic-
# composition grains and medium-composition grains, respectively.
huemulMafic <- geoEllipsoidDataFromAvizoFile("data/HU150205mafic.tsv", doNormalize=TRUE)
huemulMedium <- geoEllipsoidDataFromAvizoFile("data/HU150205medium.tsv", doNormalize=TRUE)

# Try various plots, with the mafic data set in red and the medium data set in 
# blue. It's hard to see much, but the Hsu-Nadai plot suggests that the mafic 
# data contain more near-spheres than the medium data do.
ellHsuNadaiPlot(c(huemulMafic$logA, huemulMedium$logA),
                colors=c(replicate(length(huemulMafic$logA), "red"), replicate(length(huemulMafic$logA), "blue")))
ellEqualAreaPlot(c(huemulMafic$rotation, huemulMedium$rotation),
                 c(huemulMafic$a, huemulMedium$a), 
                 colors=c(replicate(length(huemulMafic$rotation), "red"), replicate(length(huemulMafic$rotation), "blue")))
ellEqualVolumePlot(c(huemulMafic$rotation, huemulMedium$rotation), 
                   c(huemulMafic$a, huemulMedium$a), 
                   colors=c(replicate(length(huemulMafic$rotation), "red"), replicate(length(huemulMafic$rotation), "blue")),
                   simplePoints=TRUE)

# Here are the ellipsoid vector plots. Both data sets look like big clouds 
# centered at the origin. Ellipsoid data sets tend to look like that.
ellPairsPlot(huemulMafic$vector)
ellPairsPlot(huemulMedium$vector)



### TWO HYPOTHESIS TESTS ###

# In a hypothesis test, we declare a null hypothesis and an alternative 
# hypothesis. We produce a p-value, which is the probability of seeing the 
# data (or data more 'extreme' than them) if the null hypothesis is true. If 
# the p-value is less than some threshold, typically 0.05, then we reject the 
# null hypothesis.

# In this case, our null hypothesis is that the two data sets come from 
# populations with the same mean. In other words, the difference between the 
# means is zero. The alternative hypothesis is that the difference is not zero.

# We try two methods. The first method is based on a standard multivariate 
# statistics technique.
ellTwoSampleHotellingT2Inference(huemulMafic$vector, huemulMedium$vector)

# The second method is based on bootstrapping, which we have used in other 
# tutorials. This might take a minute.
huemulInf <- ellTwoSampleBootstrapInference(huemulMafic$vector, huemulMedium$vector, numBoots=10000)
huemulInf$pvalue(c(0, 0, 0, 0, 0))



### EFFECT SIZE ###

# When reporting the results of a hypothesis test, it is important to give the 
# reader enough information, so that she can understand and interpret the 
# result. Whenever possible, the report should describe the 'effect size'.

# In this example, we should give some indication of how different from zero 
# the difference in means is. The easiest way is to show plots of the 
# bootstrapped differences.
ellPairsPlot(huemulInf$bootstraps)

# A question for you to ponder: Exactly where in the plot above do we see that 
# the difference in means is not zero?

# Of course, the problem with plots of log-ellipsoid vectors is that they're 
# hard to relate to magnitude and orientation. So let's plot those aspects.
huemulInfElls <- lapply(huemulInf$bootstraps, ellEllipsoidFromVector)
ellHsuNadaiPlot(lapply(huemulInfElls, function(e) e$logA))
ellEqualAreaPlot(lapply(huemulInfElls, function(e) e$rotation), 
                 lapply(huemulInfElls, function(e) e$a))
ellEqualVolumePlot(lapply(huemulInfElls, function(e) e$rotation), 
                   lapply(huemulInfElls, function(e) e$a), simplePoints=TRUE)

# Remember that we're not actually plotting ellipsoids here. We're plotting 
# differences between ellipsoids, re-interpreted as ellipsoids. The plots 
# above give you an idea of what ellipsoid would have to be 'added' to the 
# medium-composition ellipsoids to get the mafic-composition ellipsoids. But 
# what does it mean to 'add' ellipsoids? It means convert them to log-
# ellipsoid vectors, add those, and convert back to ellipsoids. So these plots 
# of differences are, I would say, irredeemably abstract.



### CONCLUSION ###

# Two-sample inference techniques let you compare two data sets. In 
# particular, you can investigate whether their populations have the same mean.


