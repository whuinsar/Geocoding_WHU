#pragma once

#include <vector>
#include <string>

#include <boost/filesystem.hpp>
#include <gdal.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fmt/format.h>

namespace fs = boost::filesystem;

namespace WHU {
    // 0 - 99: Covering Data
    enum class DataSource {
        THU_Gong_COVER = 0,
        ESA_COVER,
        PLACEHOLDER = 100
    };
    
    struct TifEntry {
        std::string path;
        double west, east, south, north;

        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive& ar, const unsigned int) {
            ar& path
                & west& east& south& north;
        };
    };

    struct TifIndex
    {
        DataSource type;
        std::vector<TifEntry> entries;

        // TODO: Optimize this search procedure.
        std::vector<std::string> search(double west, double east, double south, double north) {
            std::vector<std::string> results;
            for (auto& entry : entries) {
                bool horizontal_intersection = !(entry.east < west || entry.west > east),
                    vertical_intersection = !(entry.north < south || entry.south > north);
                bool intersection = horizontal_intersection && vertical_intersection;
                if (intersection) {
                    results.push_back(entry.path);
                }
            };

            return results;
        };

        void add_entry(const fs::path& tif_path) {
            static bool registered = false;
            if (!registered) {
                GDALAllRegister();
            }
            
            auto dataset = (GDALDataset*)GDALOpen(tif_path.string().c_str(), GA_ReadOnly);
            int nBands = dataset->GetRasterCount();
            double transform[6]; // [0] lon, [1] lon delta, [3] lat, [5] lat delta
            if (dataset->GetGeoTransform(transform)) {
                throw std::runtime_error("Error occured when opening dataset " + tif_path.string() + "\n");
            };

            auto width = dataset->GetRasterXSize(), height = dataset->GetRasterYSize();
            TifEntry tmp{
                tif_path.filename().string(), 
                transform[0],
                transform[0] + transform[1] * width,
                transform[3] + transform[5] * height,
                transform[3],
            };

            entries.push_back(tmp);
            fmt::print("Add entry {} to list.\n", tif_path.string());
            dataset->Close();
            return;
        }

        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive& ar, const unsigned int) {
            ar& entries &type;
        };
    };
    
}