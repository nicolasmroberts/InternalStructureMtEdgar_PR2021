


# Copyright 2017 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# This file is nothing more than a convenient way to load our entire R library and its dependencies (excluding the C part).



# Load external libraries. For example, expm must be loaded before def.R, because the latter defines defExp to be expm.
library("rgl")
library("fields")
library("MASS")
library("ICSNP")
# library("FRB")
library("expm")
library("Directional")
library("pracma")

# Load our geologyGeometry library, excluding the C part.
source("geologyGeometry/library/miscellany.R")
source("geologyGeometry/library/ray.R")
source("geologyGeometry/library/rayFisher.R")
source("geologyGeometry/library/rayRegression.R")
source("geologyGeometry/library/rayPlot.R")
source("geologyGeometry/library/ray2D.R")
source("geologyGeometry/library/line.R")
source("geologyGeometry/library/lineUniform.R")
source("geologyGeometry/library/lineWatson.R")
source("geologyGeometry/library/lineBingham.R")
source("geologyGeometry/library/lineRegression.R")
source("geologyGeometry/library/lineWellner.R")
source("geologyGeometry/library/linePlot.R")
source("geologyGeometry/library/line2D.R")
source("geologyGeometry/library/rot.R")
source("geologyGeometry/library/rotRancourt.R")
source("geologyGeometry/library/rotUniform.R")
source("geologyGeometry/library/rotFisher.R")
source("geologyGeometry/library/rotRegression.R")
source("geologyGeometry/library/rotPlot.R")
source("geologyGeometry/library/ori.R")
source("geologyGeometry/library/oriRegression.R")
source("geologyGeometry/library/oriPlot.R")
source("geologyGeometry/library/ell.R")
source("geologyGeometry/library/ellSPO.R")
source("geologyGeometry/library/ellPlot.R")
source("geologyGeometry/library/def.R")
source("geologyGeometry/library/defHomogeneous.R")
source("geologyGeometry/library/defJeffery.R")
source("geologyGeometry/library/defEshelby.R")
source("geologyGeometry/library/geology.R")


