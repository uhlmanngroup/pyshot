# distutils : language = c++
# distutils: sources = src/shot_descriptor.cpp
# cython: language_level = 3
from libcpp.vector cimport vector
import numpy as np
cimport numpy as np
from libcpp cimport bool

cdef extern from "include/shot_descriptor.h":

    vector[vector[double]]  calc_shot(
                   const vector[vector[double]] vertices,
                   const vector[vector[int]] faces,
                   double radius,
                   double localRFradius,
                   int minNeighbors,
                   int bins,
                   bool double_volumes,
                   bool use_interpolation,
                   bool use_normalization,
    )

cpdef get_descriptors(
        np.ndarray[double, ndim=2] vertices,
        np.ndarray[int, ndim=2] faces,
        double radius,
        double local_rf_radius,
        int min_neighbors = 3,
        int n_bins = 20,
        bool double_volumes=True,
        bool use_interpolation=True,
        bool use_normalization=True,
):
    """
    Returns the SHOT descriptors of a point cloud 
    
    Parameters
    ------------
    vertices : (n, 3) float
      Array of vertex locations.
    faces : (m, 3)  int
      Array of triangular faces.
    radius: float
      The radius for querying neighbours.
    local_rf_radius: float
      The radius for the Reference Frame.
    min_neighbors: int
      The minimum number of neighbours to use. 
    n_bins:
      The number of bins for the histogram
    double_volumes: bool
      Double the volumes to use for the computations.
    use_interpolation: bool
      Use interpolation during computations.    
    use_normalization: bool
      Normalize during computations.    
      
    Returns
    ----------
    descr: (n, d) float
      Array containing the d SHOT descriptors for the n points,
      where d = 16 * (n_bins + 1) * (double_volumes + 1).

    """

    descr = calc_shot(vertices,
                      faces,
                      radius,
                      local_rf_radius,
                      min_neighbors,
                      n_bins,
                      double_volumes,
                      use_interpolation,
                      use_normalization)

    return np.array(descr)
