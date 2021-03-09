# distutils : language = c++
# distutils: sources = src/shot_descriptor.cpp
# cython: language_level = 3
from libcpp.vector cimport vector

cdef extern from "include/shot_descriptor.h":

    vector[vector[double]]  calc_shot(
                   vector[vector[double]] vertices,
                   vector[vector[int]] faces,
                   double radius,
                   double localRFradius,
                   int minNeighbors,
                   int bins)

cpdef get_shot_descriptor(vertices,
                         faces,
                         radius,
                         local_rf_rradius,
                         min_neighbors,
                         min):
    """
    Get shot descriptors
    
    :param vertices: 
    :param faces: 
    :param radius: 
    :param local_rf_rradius: 
    :param min_neighbors: 
    :param min: 
    :return: 
    """

    return calc_shot(vertices,
                         faces,
                         radius,
                         local_rf_rradius,
                         min_neighbors,
                         min)
