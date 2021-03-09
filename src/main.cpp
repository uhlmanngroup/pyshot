#include "shot_descriptor.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

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
