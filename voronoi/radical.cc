// Radical Voronoi tessellation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

// add the following line in Line 2084 in cell.cc
// case 'd': fprintf(fp,"%g",4.0/3.0*3.1415926*r*r*r/volume());break;

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double x_min=-100,x_max=100;
const double y_min=-102,y_max=102;
const double z_min=0,z_max=160;

// Set up the number of blocks that the container is divided
// into.
const int n_x=30,n_y=30,n_z=30;

int main() {
	// Create a container for polydisperse particles using the same
	// geometry as above. Import the polydisperse test packing and
	// output the Voronoi radical tessellation in gnuplot and POV-Ray
	// formats.
	container_poly conp(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			true,false,false,8);
	conp.import("cuboid_poly");
	conp.draw_cells_pov("cuboid_poly_v.pov");
	conp.draw_particles_pov("cuboid_poly_p.pov");

	conp.print_custom("%i %q %r %v %c %d","packing.custom");
}
