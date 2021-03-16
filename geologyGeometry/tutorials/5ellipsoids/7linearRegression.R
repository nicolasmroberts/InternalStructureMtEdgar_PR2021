


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

# In statistics, regression is a technique for explaining the variation in one 
# variable in terms of another variable. This tutorial explains scalar-
# ellipsoid  regression --- that is, fitting a curve of ellipsoids 
# parametrized by a scalar. One can use this technique to smooth a data set, 
# predict missing values, etc.



### EXPLORATION ###

# Load a file of AMS ellipsoids from the Troodos ophiolite, Cyprus (Titus et 
# al., in prep). Each ellipsoid is tagged with its geographic location, 
# measured east and north from a benchmark location. So we plot the stations 
# in map view. I have superimposed a red dot for the town of Mandria. The 
# Arakapas fault zone strikes east from Mandria, and the Solea graben strikes 
# north from Mandria (approximately). They are believed to constitute a fossil 
# ridge-transform system. 
cypAMS <- geoEllipsoidDataFromIRMFile("data/cyprusAMS.csv", doNormalize=TRUE)
plot(x=c(cypAMS$easting, 484765), y=c(cypAMS$northing, 3858516), 
     col=c(replicate(length(cypAMS$easting), "black"), "red"),
     xlab="easting (m)", ylab="northing (m)")

# Each ellipsoid is also tagged with its rock type --- mainly gabbro (green) 
# or sheeted dike complex (blue), but also a few others (black).
cypAMS$color <- "black"
cypAMS$color[grepl("gabbro", cypAMS$rock.type)] <- "green"
cypAMS$color[grepl("dike", cypAMS$rock.type)] <- "blue"
plot(x=c(cypAMS$easting, 484765), y=c(cypAMS$northing, 3858516), 
     col=c(cypAMS$color, "red"),
     xlab="easting (m)", ylab="northing (m)")

# To avoid confounding our results by rock type, we're going to focus on the 
# gabbros. To avoid dealing with different kinematics in differing parts of 
# the ridge-transform system, we're going to focus on the 'inside corner' NE 
# of Mandria. We are left with 247 ellipsoids at 26 stations.
cypInsideGabbroSelectors <- grepl("gabbro", cypAMS$rock.type) & cypAMS$easting > 484765 & cypAMS$northing > 3858516
cypInsGab <- cyp[cypInsideGabbroSelectors,]
nrow(cypInsGab)
length(levels(factor(cypInsGab$site.number)))

# Let's convert eastings and northings from m to km, and let's measure 
# relative to Mandria. The stations fill out a 20-km-by-10-km region pretty 
# well, except in the northwest part.
cypEs <- (cypInsGab$easting - 484765) / 1000
cypNs <- (cypInsGab$northing - 3858516) / 1000
plot(x=cypEs, y=cypNs, xlab="easting (km)", ylab="northing (km)")

# After playing with the data for a while, we decide to inspect the ellipsoid 
# magnitudes colored by northeasting. A vague rainbow pattern emerges: reds 
# and oranges on the left, yellows and greens near the vertex, and blues and 
# magentas on the right.
cypNEs <- (cypEs + cypNs) / sqrt(2)
ellHsuNadaiPlot(cypInsGab$logA, es=0.15, colors=hues(cypNEs))

# So maybe we should investigate whether the ellipsoids depend on 
# northeasting. One of the best places to inspect variation is a plot of the 
# first two principal components.
cypPCA <- ellPrincipalComponentAnalysis(cypInsGab$vector)
cypInsGab$pca <- lapply(cypInsGab$vector, function(v) as.numeric(t(cypPCA$rotation) %*% (v - cypPCA$center)))
ellVectorPlot(c(1, 2), cypInsGab$pca, colors=hues(cypNEs))

# That plot is not a slam dunk, but let's proceed, if only for the sake of the 
# learning.



### LINEAR REGRESSION TO FIT A POLYNOMIAL ###

# R is great in some ways, but I don't like some of its design decisions. So I 
# try to hide (what I regard as) its ugliness from you whenever possible. In 
# the case of linear regression, I haven't found any way to hide the ugliness. 
# So here we go.

# First we package the AMS log-ellipsoid vectors in a certain way. Then we 
# declare a linear model that tries to predict those vectors based on a 
# degree-4 polynomial in northeasting. (There is no physical motivation for 
# this curve. It's just mathematically tractable.)
cypVs <- t(simplify2array(cypInsGab$vector))
cypRegr <- lm(cypVs ~ 1 + cypNEs + I(cypNEs^2) + I(cypNEs^3) + I(cypNEs^4))

# Compute the predictions and convert them back into ellipsoids.
cypPreds <- predict(cypRegr, data.frame(cypNEs=seq(from=min(cypNEs), to=max(cypNEs), length.out=1000)))
cypPreds <- lapply(1:nrow(cypPreds), function(i) ellEllipsoidFromVector(cypPreds[i,]))

# The Hsu-Nadai plot displays the left-vertex-right tendency that we suspected 
# earlier, but in a complicated way.
ellHsuNadaiPlot(cypInsGab$logA, colors=hues(cypNEs), es=0.15,
                curves=list(lapply(cypPreds, function(e) e$logA)))
ellEqualVolumePlot(cypInsGab$rotation, cypInsGab$a, 
                   rotCurves=list(lapply(cypPreds, function(e) e$rotation)), 
                   aCurves=list(lapply(cypPreds, function(e) e$a)), 
                   colors=hues(cypNEs), curveWidth=3)

# For reasons that will become clear in a moment, let's view the regression's 
# predictions for the first and fifth components. And then the other three.
ellVectorPlot(c(1, 5), cypInsGab$vector, colors=hues(cypNEs),
              curves=list(lapply(cypPreds, function(e) e$vector)))
ellVectorPlot(c(2, 3, 4), cypInsGab$vector, colors=hues(cypNEs),
              curves=list(lapply(cypPreds, function(e) e$vector)), simplePoints=TRUE, curveWidth=3)

# Here's a list of the best-fit coefficients. The five columns correspond to 
# the five components of the log-ellipsoid vectors. The rows correspond to the 
# terms in the polynomial that we defined above.
coefficients(cypRegr)

# The 'summary' command gives more detailed results for each of the five 
# components (which it calls Y1, Y2, Y3, Y4, Y5). We'll discuss just a couple 
# of aspects.
summary(cypRegr)

# First, many of the coefficients are found to be statistically significantly 
# different from 0. This is especially true in the first and fifth components. 
# That's why we focused on those components in the plots above.

# Second, R^2 varies between 0.04 and 0.32. The R^2 statistic quantifies how 
# much of the variation in the data is explained by the regression. It varies 
# between 0 and 1, with 0 meaning 'no explanation' and 1 meaning 'perfect 
# explanation' (i.e. the data lying exactly on the prediction curve). What 
# constitutes a 'good' R^2 varies by discipline. In our experience, R^2 values 
# between 0.2 and 0.7 are pretty common in geology.



### USING THE RESULTS ###

# How do we use the regression, practically? First, remember that our data set 
# has a big hole in the NW part of the inside corner, for example at easting 5 
# km, northing 10 km.
plot(x=cypEs, y=cypNs, xlab="easting (km)", ylab="northing (km)")

# We can use the regression to predict what the AMS would be there.
cypPred <- as.numeric(predict(cypRegr, data.frame(cypNEs=((5 + 10) / sqrt(2)))))
cypPred <- ellEllipsoidFromVector(cypPred)
cypPred

# Second, suppose that we believe the regression in its first and fifth 
# components but not the others. We can smooth the data by replacing the first 
# and fifth components with their predictions.
cypSmoothed <- predict(cypRegr, data.frame(cypNEs=cypNEs))
cypSmoothed <- lapply(1:length(cypNEs), 
                      function(i) c(cypSmoothed[[i, 1]], cypInsGab$vector[[i]][2:4], cypSmoothed[[i, 5]]))
cypSmoothed <- lapply(cypSmoothed, ellEllipsoidFromVector)

# Here are various plots of the smoothed data. The (1, 5) vector plot looks 
# really weird, because those are exactly the components that we smoothed. In 
# other words, that plot is just showing the prediction curve.
ellHsuNadaiPlot(lapply(cypSmoothed, function(e) e$logA), colors=hues(cypNEs), es=0.15)
ellEqualVolumePlot(lapply(cypSmoothed, function(e) e$rotation), lapply(cypSmoothed, function(e) e$a), colors=hues(cypNEs))
ellVectorPlot(c(1, 5), lapply(cypSmoothed, function(e) e$vector), colors=hues(cypNEs))
ellVectorPlot(c(2, 3, 4), lapply(cypSmoothed, function(e) e$vector), colors=hues(cypNEs))

# Third, we now have evidence that the ellipsoids are not 'independent and 
# identically distributed' (IID), because they depend on geography. Many 
# inference techniques assume that the data are IID. So we should not use 
# those techniques on this data set. That's good to know.

# Finally, we'd like to be able to verbally summarize the result of the 
# regression --- something along the lines of, 'for every km of northeasting, 
# such-and-such happens'. But this is difficult, because we fitted a non-
# linear curve to the data. If we'd fitted a linear curve, it would have a 
# simple 'for every km of northeasting...' meaning. But then a linear curve 
# would never have fit the data very well.



### CONCLUSION ###

# When ellipsoids are rendered into their log-ellipsoid vectors, you can use 
# standard multivariate regression techniques on them, to (try to) pick out 
# patterns among them.


