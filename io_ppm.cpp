#include <iostream>
#include <fstream>
#include "io_ppm.hpp"

using namespace std;

namespace raster_stats {


void write_ppm(const landscape_t& raster, const char* filename) {
	ofstream out(filename);
	// Write the width first, then the height.
	out << "P5 " << raster.size1() << " " << raster.size2() << " 16" << endl;
	for (size_t row_idx=0; row_idx<raster.size1(); row_idx++ ) {
		size_t rev_row_idx = raster.size2()-row_idx-1;
		for (size_t col_idx=0; col_idx<raster.size2(); col_idx++) {
			//cout << rev_row_idx << " " << col_idx << endl;
			out << (int) raster(rev_row_idx,col_idx) << " ";
		}
		out << endl;
	}
}


}
