#include <math.h>
#include <vector>
#include "shot_descriptor.h"

#include <fstream>
#include <iostream>

int main(int argc, char ** argv) {
	// Create mesh structure

	mesh_t mesh;
    int nv = 1;
    std::vector< vec3d<double> > V(nv);

        vec3d<double>* b1 = new vec3d<double> (0.0, 0.0, 0.0);
//        vec3d<double> b2();//0.0, 1.0, 0.0);
//        vec3d<double> b3();//0.0, 0.0, 1.0);
//        vec3d<double> b4();//1.0, 0.0, 0.0);
//
        V.push_back(*b1);
//        V.push_back(b2);
//        V.push_back(b3);
//        V.push_back(b4);
//
//		mesh.put_vertices(V);
//        mesh.add_triangle(0, 1, 2);
//        mesh.add_triangle(0, 2, 3);
//        mesh.add_triangle(0, 1, 3);
//        mesh.add_triangle(1, 2, 3);

//	mesh.calc_normals();
//
	// Compute SHOT descriptors at the desired point indices

	unibo::SHOTParams params;
	params.radius = 12;
	params.localRFradius = params.radius;
	params.minNeighbors = 4;
	params.bins = 20;

	unibo::SHOTDescriptor sd(params);
	const size_t sz = sd.getDescriptorLength();

    int np = 42;

	std::cout << "Computing SHOTs on " << np << " points... " << std::flush;

	for (size_t i=0; i<np; ++i)
	{
		unibo::shot s;
		//sd.describe(mesh, static_cast<int>(idx[i]-1), s);
		//for (size_t j=0; j<sz; ++j)
		//s	D[i*sz+j] = s(j);
	}

	std::cout << "done." << std::endl;
}
