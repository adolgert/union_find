#ifndef _IO_HDF_HPP_
#define _IO_HDF_HPP_ 1


#include <string>
#include <sstream>
#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/filesystem.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/uuid/name_generator.hpp>
#include <boost/uuid/string_generator.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include "H5Cpp.h"
#include "H5File.h"
#include "H5Group.h"
#include "H5FaccProp.h"
#include "hdf5.h"

// Mac address
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <net/if.h>


// versions
#include <boost/version.hpp>
#include "tbb/tbb_stddef.h"
#include "python2.7/patchlevel.h"
#include "raster_version.hpp"
#include "tiffvers.h"


/*! Write timing data.
 *  That's an array of number of iterations and total time, for each repitition.
 *  The metadata is: compilation flags, command line arguments, version
 *  of libraries, svn info for main code, maybe some indication of what
 *  is timed.
 */


namespace raster_stats
{
    class timing_file
    {
        std::string filename_;
        H5::H5File file_;
        boost::uuids::uuid timing_id_;
        boost::uuids::uuid build_id_;
        boost::uuids::uuid machine_id_;
        bool is_ok_;
        boost::uuids::uuid namespace_uuid_;
    public:
        /*! Create a file for timing.
         *  After the file is created, check is_ok() to see if there were errors.
         */
        timing_file(const std::string& filename,
                    const boost::program_options::basic_parsed_options<char>& parsed_options)
            : filename_(filename)
        {
            boost::uuids::string_generator sgen;
            namespace_uuid_=sgen("6ba7b810-9dad-11d1-80b4-00c04fd430c8");

            is_ok_=false;

            // The HDF5 C++ interface uses exceptions to tell you when you can and cannot
            // open a group. That means we have to use exceptions two ways. The large try-catch
            // is designed to find when the file succeeded or failed to create, but smaller
            // try-catch blocks within that are used as part of the control flow.

            try {
                boost::filesystem::path filepath(filename);
                if (boost::filesystem::exists(filepath)) {
                    // You can't use H5F_ACC_RDWR flag. Must instead call openFile
                    // whose default is RDWR. Not documented.
                    //file_.openFile(filename,H5F_ACC_RDWR|H5F_ACC_DEBUG); // 0x1u|H5F_ACC_DEBUG);
                    // And openFile doesn\'t do what it should, so use H5Fopen.
                    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                    file_.setId(file_id);
                } else {
                    file_=H5::H5File(filename,H5F_ACC_TRUNC|H5F_ACC_DEBUG);
                }

                write_machine();

                write_build();


                H5::Group timing_group;
                try {
                    timing_group = H5::Group(file_.openGroup("/timing"));
                    // catch Exception b/c missing group may through Group or File exceptions.
                } catch (H5::Exception error) {
                    timing_group = H5::Group(file_.createGroup("/timing"));
                }

                // Name the dataset with its uuid.
                std::stringstream dataset_base;
                dataset_base << "/timing/" << timing_id_;
                H5::Group instance_group(file_.createGroup(dataset_base.str().c_str()));

                // In the timing run, store the machine and build.
                std::stringstream build_dir;
                build_dir << "/build/" << build_id_;
                write_attribute(instance_group,"build",build_dir.str());

                std::stringstream machine_dir;
                machine_dir << "/machine/" << machine_id_;
                write_attribute(instance_group,"machine",machine_dir.str());

                // When did this run start?
                boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
                std::stringstream whenstr;
                whenstr << now;
                write_attribute(instance_group,"start time",whenstr.str());

                // Save command-line options for this run.
                auto pidx=parsed_options.options.begin();
                while (pidx!=parsed_options.options.end()) {
                    std::string opt_name(pidx->string_key);
                    if (pidx->position_key>=0) {
                        std::stringstream pos_str;
                        pos_str << "pos" << pidx->position_key;
                        opt_name=pos_str.str();
                    }
                    std::stringstream val_str;
                    for (size_t vidx=0; vidx<pidx->value.size(); vidx++) {
                        val_str << pidx->value[vidx] << std::endl;
                    }
                    write_attribute(instance_group,opt_name,val_str.str());
                    pidx++;
                }

                is_ok_=true;
            }
            catch (H5::FileIException error) {
                std::cerr << "FileIException" << std::endl;
                error.printError();
                return;
            }
            // catch failure caused by the DataSpace operations
            catch( H5::DataTypeIException error )
            {
               std::cerr << "DataTypeIException" << std::endl;
               error.printError();
               return;
            }

        }

        bool is_ok() { return is_ok_; }


        void store_test(const std::vector<boost::array<size_t,2>>& data,
                        const std::string& timing_name)
        {
            if (!is_ok()) {
                std::cerr << "Skipping writing " << timing_name << " because file isn\'t open."
                        << std::endl;
                return;
            }
            std::cout << "Writing " << timing_name << " to " << filename_ << std::endl;

            try {

                hsize_t dimsf[2];
                dimsf[0] = data.size();
                dimsf[1] = 2;
                int rank = 2;
                H5::DataSpace data_space( rank, dimsf );

                H5::IntType data_type( H5::PredType::NATIVE_UINT64 );
                data_type.setOrder( H5T_ORDER_LE );

                std::stringstream dataset_name;
                dataset_name << "/timing/" << timing_id_ << "/" << timing_name;
                H5std_string DATASET_NAME( H5std_string(dataset_name.str()) );

                H5::DataSet dataset = file_.createDataSet( DATASET_NAME,
                                                                data_type,
                                                                data_space );

                dataset.write( &data[0], H5::PredType::NATIVE_UINT64 );
            }

            catch (H5::FileIException error) {
                std::cerr << "FileIException" << std::endl;
                error.printError();
                return;
            }
           // catch failure caused by the DataSet operations
           catch( H5::DataSetIException error )
           {
              std::cerr << "DataSetIException" << std::endl;
              error.printError();
              return;
           }

           // catch failure caused by the DataSpace operations
           catch( H5::DataSpaceIException error )
           {
              std::cerr << "DataSpaceIException" << std::endl;
              error.printError();
              return;
           }
    
           // catch failure caused by the DataSpace operations
           catch( H5::DataTypeIException error )
           {
              std::cerr << "DataTypeIException" << std::endl;
              error.printError();
              return;
           }
        }



        void write_machine()
        {
            struct ifreq buffer;
            // b8:ac:6f:a2:5e:f6
            std::stringstream ethstr;

			#ifdef _LINUX
            int fd = socket(PF_INET,SOCK_DGRAM,0);
            memset(&buffer,0x00,sizeof(buffer));
            strcpy(buffer.ifr_name,"eth0");
            if (fd>=0 && 0==ioctl(fd,SIOCGIFHWADDR, &buffer)) {
                ethstr << std::hex << std::setw(2) << std::setfill('0')
                    << static_cast<int>(static_cast<unsigned char>(buffer.ifr_hwaddr.sa_data[0]));

                for (size_t si=1; si<6; si++) {
                    char v = buffer.ifr_hwaddr.sa_data[si];
                    ethstr << ":" //<< std::hex << std::setw(2) << std::setfill('0')
                        << static_cast<int>(static_cast<unsigned char>(v));
                }
            } else {
                ethstr << "de:ad:de:ad:be:ef";
            }
            close(fd);
			#endif
            ethstr << "and more machine details to increase hash length";


            boost::uuids::name_generator gen(namespace_uuid_);
            machine_id_=boost::uuids::uuid(gen(ethstr.str()));

            H5::Group machine_group;
            try {
                std::cout << "looking for machine group" << std::endl;
                machine_group = file_.openGroup("/machine");
                std::cout << "found already existing machine group" << std::endl;
            }
            catch( H5::Exception error ) {
                machine_group = file_.createGroup("/machine");
                std::cout << "created new /machine group" << std::endl;
            }

            std::stringstream machine_base;
            machine_base << "/machine/" << machine_id_;

            try {
                std::cout << "machine name:" << machine_base.str() << std::endl;
                H5::Group host = file_.openGroup(machine_base.str());
                std::cout << "found the machine already there" << std::endl;
            }
            catch( H5::Exception error ) {
                std::cout << "about to create a new machine group for this machine" << std::endl;
                H5::Group instance_group=H5::Group(file_.createGroup(machine_base.str().c_str()));
                std::cout << "created a new machine group" << std::endl;;

                write_attribute(instance_group,"ethernet mac address",ethstr.str());
            }
        }



       void write_build()
        {
           std::stringstream version_string;
           version_string << BOOST_VERSION << PY_MAJOR_VERSION << PY_MINOR_VERSION
                   << PY_MICRO_VERSION << TIFFLIB_VERSION << TBB_VERSION_MAJOR
                   << TBB_VERSION_MINOR
#ifdef __GNUC__
                   << __GNUC__ << __GNUC_MINOR__ << __GNUC_PATCHLEVEL__
#endif
                   << RASTER_STATS_VERSION << RASTER_STATS_CFG
                   << RASTER_STATS_COMPILE_TIME;


           boost::uuids::name_generator gen(namespace_uuid_);
           build_id_=boost::uuids::uuid(gen(version_string.str()));

            try {
                // Make a /build if it does not exist.
                try {
                    auto build_group = H5::Group(file_.openGroup("/build"));
                    std::cout << "found the build group" << std::endl;
                }
                catch (H5::Exception error) {
                    auto build_group = H5::Group(file_.createGroup("/build"));
                    std::cout << "created a build group" << std::endl;
                }

                // Name the dataset with its uuid.
                std::stringstream dataset_base;
                dataset_base << "/build/" << build_id_;

                H5::Group instance_group;
                // If this particular build exists, no cause to duplicate it.
                try {
                    instance_group = file_.openGroup(dataset_base.str());
                    std::cout << "found this particular build" << std::endl;
                }
                catch (H5::Exception error) {
                    instance_group = H5::Group(file_.createGroup(dataset_base.str()));
                    std::cout << "created a new build" << std::endl;

                    std::stringstream bstr;
                    bstr << BOOST_VERSION;
                    write_attribute(instance_group,"Boost version",bstr.str());

                    std::stringstream pystr;
                    pystr << PY_MAJOR_VERSION << "." << PY_MINOR_VERSION << "." << PY_MICRO_VERSION;
                    write_attribute(instance_group,"Python version",pystr.str());

                    std::stringstream tbbstr;
                    tbbstr << TBB_VERSION_MAJOR << "." << TBB_VERSION_MINOR;
                    write_attribute(instance_group,"TBB version",tbbstr.str());

    #ifdef __GNUC__
                    std::stringstream gnustr;
                    gnustr << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__;
                    write_attribute(instance_group,"GNUC version",gnustr.str());
                    write_attribute(instance_group,"compiler","GNUC version");
    #endif

                    write_attribute(instance_group,"Raster stats version",
                            RASTER_STATS_VERSION);
                    write_attribute(instance_group,"Raster stats configuration",
                            RASTER_STATS_CFG);
                    write_attribute(instance_group,"Raster stats compile time",
                            RASTER_STATS_COMPILE_TIME);
                }
            }
            catch (H5::FileIException error) {
                std::cerr << "FileIException" << std::endl;
                error.printError();
                return;
            }
            // catch failure caused by the DataSpace operations
            catch( H5::DataTypeIException error )
            {
               std::cerr << "DataTypeIException" << std::endl;
               error.printError();
               return;
            }
        }


       void write_attribute(H5::Group& group, std::string name, std::string value)
       {
           H5::DataSpace attr_dataspace(H5S_SCALAR);
           H5::StrType base_t(H5::PredType::C_S1, H5T_VARIABLE);
           const H5std_string strwritebuf(value);

           auto bver=group.createAttribute(name,base_t,attr_dataspace);
           bver.write(base_t, value);

       }




       std::string uuid_group(boost::uuids::uuid u) {
           std::string s = to_string(u);
           size_t pos = std::string::npos;
           size_t found;
           while ((found=s.find('-',pos)) != std::string::npos) {
               s.erase(found,1);
           }
           return s;
       }
    };

}


#endif // _IO_HDF_HPP_
