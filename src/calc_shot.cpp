#include "shot_descriptor.h"
#include <math.h>
#include <vector>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

int main(int argc, char **argv) {
  mesh_t mesh;

  // Parse a simple OFF file without any header
  std::ifstream infile(argv[1]);
  std::cout << "Reading from " << argv[1] << std::endl;

  int nv, nf, ne;
  infile >> nv >> nf >> ne;

  std::vector<vec3d<double>> V(nv);
  {
    double x, y, z;
    for (int i = 0; i < nv; i++) {
      infile >> x >> y >> z;
      vec3d<double> b(x, y, z);
      V[i] = b;
    }
  }
  mesh.put_vertices(V);
  {
    int a, b, c;
    for (int i = 0; i < nf; i++) {
      infile >> a >> b >> c;
      mesh.add_triangle(a, b, c);
    }
  }

  mesh.calc_normals();

  // Compute SHOT descriptors at the desired point indices
  unibo::SHOTParams params;
  params.radius = 10000;
  params.localRFradius = params.radius;
  params.minNeighbors = 4;
  params.bins = 20;

  unibo::SHOTDescriptor sd(params);
  const size_t sz = sd.getDescriptorLength();

  int np = nv;

  for (int i = 0; i < nv; i++) {
    unibo::shot s;
    sd.describe(mesh, i, s);
    // for (size_t j=0; j<sz; ++j)
    //     std::cout << "i=" << i << " j=" << j << " " << s(j) << std::endl;
  }

  std::cout << "Computed " << nv * sz << " for SHOT descriptors" << std::endl;
}
