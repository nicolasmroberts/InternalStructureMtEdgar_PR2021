


# Copyright 2016 Joshua R. Davis
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.



# A line is expressed as a unit 3D vector in Cartesian coordinates, as an R vector u = c(x, y, z) where x^2 + y^2 + z^2 == 1. It is crucial to remember that u and -u represent the same line. Many of the line functions here are implemented in terms of the ray functions of rays.R. In particular, when lines are tightly concentrated, you can often ignore their 'negative copies' and treat them as rays, with no appreciable effect on the statistics.



### UNIFORM DISTRIBUTION ###

#' Uniformly random lines.
#' 
#' @param n A real number (positive integer) or NULL.
#' @return If n is NULL, then a single line. If n is a positive integer, then a list of n lines.
lineUniform <- function(n=NULL) {
  if (is.null(n))
    lower(rayUniform())
  else
    lapply(rayUniform(n), lower)
}

# Bingham test of uniformity, Mardia and Jupp (2000, p. 232)
# Gine test of uniformity, p. 233


