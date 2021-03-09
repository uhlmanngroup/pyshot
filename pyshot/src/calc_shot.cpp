#include "shot_descriptor.h"
#include <math.h>
#include <vector>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

std::vector<std::vector<double > >  calc_shot(
               std::vector<std::vector<double> > vertices,
               std::vector<std::vector<int> > faces,
               double radius,
               double localRFradius,
               int minNeighbors,
               int bins
)
{
  mesh_t mesh;
  int nv = vertices.size();
  int nf = faces.size();

  std::vector<vec3d<double>> V(nv);

  for (int i = 0; i < nv; i++) {
    vec3d<double> b(vertices[i][0], vertices[i][1], vertices[i][2]);
    V[i] = b;
  }
  mesh.put_vertices(V);
  for (int i = 0; i < nf; i++) {
    mesh.add_triangle(faces[i][0], faces[i][1], faces[i][2]);
  }

  mesh.calc_normals();

  unibo::SHOTParams params;
  params.radius = radius;
  params.localRFradius = localRFradius;
  params.minNeighbors = minNeighbors;
  params.bins = bins;

  unibo::SHOTDescriptor sd(params);
  const size_t sz = sd.getDescriptorLength();

  std::vector<std::vector<double > > descriptors(nv, std::vector<double>(sz));

  for (int i = 0; i < nv; i++) {
    unibo::shot s;
    sd.describe(mesh, i, s);
    for (size_t j=0; j < sz; j++) {
        descriptors[i][j] = (double) s(j);
    }
  }
  return descriptors;
}

int main(int argc, char **argv) {

  // Parse a simple OFF file without any header
  std::ifstream infile(argv[1]);
  std::cout << "Reading from " << argv[1] << std::endl;

  int nv, nf, ne;
  infile >> nv >> nf >> ne;

  std::vector<std::vector<double > > vertices(nv, std::vector<double>(3));
  std::vector<std::vector<int > > faces(nf, std::vector<int>(3));

  double x, y, z;
  for (int i = 0; i < nv; i++) {
    infile >> vertices[i][0] >> vertices[i][1] >> vertices[i][2];
  }

  int a, b, c;
  for (int i = 0; i < nf; i++) {
    infile >> faces[i][0] >> faces[i][1] >> faces[i][2];
  }

  std::vector<std::vector<double > > descriptors = calc_shot(vertices, faces,
                                                    10000, 10000, 4, 20);

  size_t sz = sizeof(descriptors) / nv;
//  const int n_descriptors = descriptors[0].size();
//  for (int i = 0; i < nv; i++) {
//    for (size_t j=0; j < n_descriptors; j++) {
//        std::cout << descriptors[i][j] << std::endl;
//    }
//  }
  std::cout << "Done" << std::endl;
}
