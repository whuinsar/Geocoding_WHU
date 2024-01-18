#pragma once

#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace WHU {
    // This interface only works under the assumption that &ISerializable = &T,
    // it should be the first entry in inheritance list, and would be better if 
    // T contains no virtual functions.
    template <typename T>
    struct ISerializable {
        static T fromBoostBinary(const boost::filesystem::path& input_filename) {
            T tmp;
            std::ifstream input_file(input_filename.string(), std::ios_base::binary);
            boost::archive::binary_iarchive ia(input_file);
            ia >> tmp;
            input_file.close();
            return tmp;
        };

        static T fromBoostText(const boost::filesystem::path& input_filename) {
            T tmp;
            std::ifstream input_file(input_filename.string());
            boost::archive::text_iarchive ia(input_file);
            ia >> tmp;
            input_file.close();
            return tmp;
        };
        void toBoostBinary(const boost::filesystem::path& output_filename) {
            std::ofstream output_file(output_filename.string(), std::ios_base::binary);
            boost::archive::binary_oarchive oa(output_file);
            // A bizzard workaround.
            T* tmp = reinterpret_cast<T*>(this);
            oa << (*tmp);
            output_file.close();
        };
        
        void toBoostText(const boost::filesystem::path& output_filename) {
            std::ofstream output_file(output_filename.string());
            boost::archive::text_oarchive oa(output_file);
            T* tmp = reinterpret_cast<T*>(this);
            oa << (*tmp);
            output_file.close();
        };

    };
}