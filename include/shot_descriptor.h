#ifndef SHOT_DESCRIPTOR_H
#define SHOT_DESCRIPTOR_H

#include "mesh.h"
#include <ostream>
#include <stdexcept>
#include <string>
#include <math.h>
#include <vector>

class invalid_mesh_descriptor : public std::logic_error {
public:
  explicit invalid_mesh_descriptor()
      : std::logic_error("Exception invalid_mesh_descriptor caught.") {}
  invalid_mesh_descriptor(const std::string &msg)
      : std::logic_error("Exception invalid_mesh_descriptor caught: " + msg) {}
};

/**
 * Generic pointwise local shape descriptor for meshes.
 */
template <typename T> class mesh_descriptor {
protected:
  T point;

public:
  explicit mesh_descriptor() {}
  mesh_descriptor(const T &p) : point(p) {}
};

/**
 * Generic pointwise indexed vector descriptor.
 */
template <typename T> struct vec_descriptor : public mesh_descriptor<int> {
protected:
  std::vector<T> d;

public:
  typedef T value_type;

  explicit vec_descriptor() : mesh_descriptor<int>() {}
  explicit vec_descriptor(size_t n) : mesh_descriptor<int>() {
    d.resize(n, T(0));
  }

  const T &operator()(size_t k) const { return d[k]; }
  T &operator()(size_t k) { return d[k]; }

  const T &at(size_t k) const { return d.at(k); }
  T &at(size_t k) { return d.at(k); }

  size_t size() const { return d.size(); }
  bool empty() const { return d.empty(); }
  void clear() { d.clear(); }
  void resize(size_t s) { d.resize(s); }
  void resize(size_t s, T v) { d.resize(s, v); }

  bool operator==(const vec_descriptor<T> &other) const {
    if (this == &other)
      return true;
    return d == other.d;
  }

  bool operator!=(const vec_descriptor<T> &other) const {
    if (this == &other)
      return false;
    return !(*this == other);
  }

  vec_descriptor<T> operator+(const vec_descriptor<T> &other) const {
    const size_t s = size();
    if (s != other.size())
      throw invalid_mesh_descriptor("Cannot sum different length descriptors.");

    vec_descriptor<T> sum;
    for (size_t k = 0; k < s; ++k)
      sum[k] = this->d[k] + other.d[k];

    return sum;
  }

  vec_descriptor<T> operator-(const vec_descriptor<T> &other) const {
    const size_t s = size();
    if (s != other.size())
      throw invalid_mesh_descriptor("Cannot sum different length descriptors.");

    vec_descriptor<T> sum;
    for (size_t k = 0; k < s; ++k)
      sum[k] = this->d[k] - other.d[k];

    return sum;
  }

  vec_descriptor<T> &operator+=(const vec_descriptor<T> &other) {
    const size_t s = size();
    if (s != other.size())
      throw invalid_mesh_descriptor("Cannot sum different length descriptors.");

    for (size_t k = 0; k < s; ++k)
      this->d[k] += other.d[k];

    return *this;
  }

  vec_descriptor<T> &operator-=(const vec_descriptor<T> &other) {
    const size_t s = size();
    if (s != other.size())
      throw invalid_mesh_descriptor("Cannot sum different length descriptors.");

    for (size_t k = 0; k < s; ++k)
      this->d[k] -= other.d[k];

    return *this;
  }

  template <typename S> vec_descriptor<T> &operator/=(const S &val) {
    const size_t s = size();
    for (size_t k = 0; k < s; ++k)
      this->d[k] /= T(val);
    return *this;
  }

  friend std::ostream &operator<<(std::ostream &s, const vec_descriptor &d) {
    s << d(0);
    for (size_t k = 1; k < d.size(); ++k) {
      s << " " << d(k);
    }
    return s;
  }
};

namespace unibo {

class SHOTDescriptor : public vec_descriptor<float> {
  typedef vec_descriptor<float> super;

public:
  typedef super::value_type value_type;

  explicit SHOTDescriptor() : super(), radius(0.0f) {}
  explicit SHOTDescriptor(size_t n) : super(n), radius(0.0f) {}

  float radius;
};

class SHOTComputer {
public:
  explicit SHOTComputer(
    float radius = 15,
    float localRFradius = 15,
    int minNeighbors = 10,
    float bins = 10,
    double doubleVolumes=true,
    double useInterpolation=true,
    double useNormalization=true
  ):
    m_radius(radius),
    m_localRFradius(radius),
    m_bins(bins),
    m_minNeighbors(minNeighbors),
    m_doubleVolumes(doubleVolumes),
    m_useInterpolation(useInterpolation),
    m_useNormalization(useNormalization),
    m_onion_husks(doubleVolumes? 32: 16)
  {}

  void describe(mesh_t &data, int feat_index, SHOTDescriptor &desc) const;

  inline int getDescriptorLength() const { return m_onion_husks * (m_bins + 1); }

private:
  // Radius of the sphere that defines the local neighborhood
  const float m_radius;

  // Radius of the support for the estimation of the local RF
  const float m_localRFradius;

  // Quantization bins for the cosine
  const int m_bins;
  const int m_minNeighbors;

  const bool m_doubleVolumes;
  const bool m_useInterpolation;
  const bool m_useNormalization;

  // Number of onion husks
  const int m_onion_husks;

/**
 * Compute the local reference frame (RF) for a given point.
 *
 * Neighbours and their associated (squared) distances are computed here.
 *
 * @param data Mesh data-structure of the point cloud.
 * @param p Index of the point.
 * @param radius Radius of the neighbourhood to consider for the RF.
 * @param X First (non normalised) vector of the reference frame
 * (modified in place)
 * @param Y Second (non normalised) vector of the reference frame
 *  (modified in place)
 * @param Z Third (non normalised) vector of the reference frame
 *  (modified in place)
 */
  void getSHOTLocalRF(mesh_t &data, int p, double radius, vec3d<double> &X,
                      vec3d<double> &Y, vec3d<double> &Z) const;

/**
 * Compute the local reference frame (RF) for a given point.
 *
 * @param data Mesh data-structure of the point cloud.
 * @param p Index of the point.
 * @param pts Vector of the neighbours' indices.
 * @param dists Vector of the squared distances associated to neighbours.
 * @param radius Radius of the neighbourhood to consider for the RF.
 * @param X First normalised vector of the reference frame (modified in place)
 * @param Y Second normalised vector of the reference frame (modified in place)
 * @param Z Third normalised vector of the reference frame (modified in place)
 */
  void getSHOTLocalRF(mesh_t &data, int p, const std::vector<int> &pts,
                      const std::vector<double> &dists, double radius,
                      vec3d<double> &X, vec3d<double> &Y,
                      vec3d<double> &Z) const;
};

} // namespace unibo

/**
 * Returns the SHOT descriptors of a mesh point cloud.
 *
 * @param vertices (n, 3) array of vertex locations.
 * @param faces (m, 3) array of triangular faces indices.
 * @param radius Radius for querying neighbours.
 * @param localRFradius Radius of the Reference Frame neighbourhood.
 * @param minNeighbors The minimum number of neighbours to use.
 * @param bins Number of bins for the histogram
 * @param double_volumes Double the maximum number of volume angular sectors
   for descriptor.
 * @param use_interpolation Use interpolation during computations.
 * @param use_normalization Normalize during computations.
 * @return (n, d) array containing the d SHOT descriptors for the n points,
   where d = 16 * (bins + 1) * (double_volumes + 1).
 */
std::vector<std::vector<double > >  calc_shot(
               const std::vector<std::vector<double> >& vertices,
               const std::vector<std::vector<int> >& faces,
               double radius,
               double localRFradius,
               int minNeighbors,
               int bins,
               bool doubleVolumes=true,
               bool useInterpolation=true,
               bool useNormalization=true
);

#endif // SHOT_DESCRIPTOR_H
