#include <map>
#include <boost/test/included/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/version.hpp>
#include "array_store.hpp"
//#include "single.hpp"
//#include "grid2d.hpp"

using namespace std;
using namespace boost::unit_test;
using namespace raster_stats;



int main(int argc, char* argv[])
{
	transform_ij tij(512);
	transform_map<transform_ij,unsigned char> data(tij,512*256);
}
