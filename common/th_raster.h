#pragma once

#include <vector>
#include <tuple>
#include <utility>

#include "th_tif_index.h"
#include "th_serialize.h"

#include "tiffio.h"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

// #define DEBUG

namespace WHU {
    template<typename T>
    struct GenericRaster : public ISerializable<GenericRaster<T>> {
        double lon_start, lat_start;
        double lon_delta, lat_delta;

        size_t width, height;

        using value_type = T;
        std::vector<T> data;

        // y for height, x, for width. (or i,j)
        T at(size_t y, size_t x) const {
            if (y < 0 || y >= height || x < 0 || x >= width) {
                if constexpr (std::is_floating_point_v<T>) {
                    return std::nan("");
                }
                else {
                    throw std::runtime_error("Raster index out of scope.\n");
                };
            }

            return data[y * width + x];
        }

        // y for height, x, for width. (or i,j)
        T& at(size_t y, size_t x) {
            if (y < 0 || y >= height || x < 0 || x >= width) {
                // In non-const version, not return anything.
                throw std::runtime_error("Raster index out of scope.\n");
            }

            return data[y * width + x];
        }

        std::tuple<double, double> coordinate(size_t y, size_t x) {
            return { lat_start + y * lat_delta, lon_start + x * lon_delta };
        };

        double distance(size_t y, size_t x, double lat, double lon) {
            auto [lat_, lon_] = coordinate(y, x);
            return std::sqrt((lat - lat_) * (lat - lat_) + (lon - lon_) * (lon - lon_));
        };

        // TODO: extract the public part between this and DEM::writeAsTiff
        void writeAsTiff(const fs::path& path) {
            if (fs::is_regular_file(path)) {
                fs::remove(path);
            }

            TIFF* tif = TIFFOpen(path.string().c_str(), "w");
            if (!tif) {
                throw std::runtime_error(
                    fmt::format("Failed to open image file: %s\n", path.string().c_str())
                );
            }

            uint16 bits_per_sample = sizeof(T) * 8;
            uint16 samples_per_pixel = 1;
            uint16 compression = COMPRESSION_NONE;
            uint16 photo_interpretation = PHOTOMETRIC_MINISBLACK;

            if (TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, width) != 1) {
                throw std::runtime_error("Failed to set image width\n");
            }

            if (TIFFSetField(tif, TIFFTAG_IMAGELENGTH, height) != 1) {
                throw std::runtime_error("Failed to set image height\n");
            }

            if (TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bits_per_sample) != 1) {
                throw std::runtime_error("Failed to set image bits_per_sample\n");
            }

            if (TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, samples_per_pixel) != 1) {
                throw std::runtime_error("Failed to set image samples_per_pixel\n");
            }

            if (TIFFSetField(tif, TIFFTAG_COMPRESSION, compression) != 1) {
                throw std::runtime_error("Failed to set image compression\n");
            }

            if (TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, photo_interpretation) != 1) {
                throw std::runtime_error("Failed to set image photo_interpretation\n");
            }

            for (size_t i = 0; i < height; i++) {
                if (TIFFWriteScanline(tif, &(data[i * width]), i) != 1) {
                    throw std::runtime_error(
                        fmt::format("Failed to write image line {}\n", i)
                    );
                }
            }
            TIFFClose(tif);

            return;
        };

        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive& ar, const unsigned int) {
            ar& lon_start& lat_start& lon_delta& lat_delta
                & width& height& data;
        };
    };

    enum class InterpolationMethod {
        Nearest,
        Bilinear
    };

    template <typename T>
    struct RasterMergeSession {

        // Pixels of output image, whose value cannot be determined from a single data source.
        struct Ambiguity {
            using base_type = std::tuple<double, double, T>; // lat, lon, val
            std::vector<base_type> candidates;
        };

        GenericRaster<T>& target;
        std::vector<Ambiguity*> ambiguities;
        std::vector<bool> mask;
        size_t number_of_pixels_without_value;
        size_t number_of_amguities;
        InterpolationMethod method;

        RasterMergeSession(GenericRaster<T>& t, InterpolationMethod m = InterpolationMethod::Nearest) :
            target(t),
            ambiguities(target.data.size(), nullptr),
            mask(target.data.size(), 0),
            number_of_pixels_without_value(target.data.size()),
            number_of_amguities(0),
            method(m)
        {
        };

        // May modify the original data.
        void from(GenericRaster<T>& source) {
            switch (method) {
            case InterpolationMethod::Nearest:
            case InterpolationMethod::Bilinear:
                fourPointsMerge(source);
                break;
            default:
                break;
            }
        };

        void fourPointsMerge(GenericRaster<T>& source) {
            assert(target.lon_delta > 0 && target.lat_delta > 0);

            auto exchange = [](double& origin, double& delta, size_t count) {
                double end = origin + (count - 1) * delta;
                std::swap(origin, end);
                delta = -delta;
            };
            if (source.lat_delta < 0) {
                for (size_t i = 0; i < source.height / 2; ++i) {
                    for (size_t j = 0; j < source.width; ++j) {
                        std::swap(source.data[i * source.width + j], source.data[(source.height - 1 - i) * source.width + j]);
                    }
                }
                exchange(source.lat_start, source.lat_delta, source.height);
                fmt::print("Vertical flip performed. lat_spacing: {} => {}.\n", -source.lat_delta, source.lat_delta);
            };
            if (source.lon_delta < 0) {
                for (size_t i = 0; i < source.height; ++i) {
                    for (size_t j = 0; j < source.width / 2; ++j) {
                        std::swap(source.data[i * source.width + j], source.data[source.height * source.width + source.width - j - 1]);
                    }
                }
                exchange(source.lon_start, source.lon_delta, source.width);
                fmt::print("Horizontal flip performed. lon_spacing: {} => {}.\n", -source.lon_delta, source.lon_delta);
            };

            int bottom_out = std::floor((source.lat_start - target.lat_start) / target.lat_delta),
                left_out = std::floor((source.lon_start - target.lon_start) / target.lon_delta),
                upper_out = std::ceil((source.lat_start + (source.height - 1) * source.lat_delta - target.lat_start) / target.lat_delta),
                right_out = std::ceil((source.lon_start + (source.width - 1) * source.lon_delta - target.lon_start) / target.lon_delta);

            int bottom_in = std::max<int>(0, bottom_out + 1),
                left_in = std::max<int>(0, left_out + 1),
                upper_in = std::min<int>(target.height - 1, upper_out - 1),
                right_in = std::min<int>(target.width - 1, right_out - 1);

            fmt::print("Inner extent : {}, {}, {}, {}\n", bottom_in, upper_in, left_in, right_in);
            for (int i = bottom_in; i <= upper_in; ++i) {
                int i_source_bottom = std::floor((target.lat_start + target.lat_delta * i - source.lat_start) / source.lat_delta),
                    i_source_upper = i_source_bottom + 1;
                for (int j = left_in; j <= right_in; ++j) {
                    if (mask[target.width * i + j] == 1) {
                        continue;
                    };

                    int j_source_left = std::floor((target.lon_start + target.lon_delta * j - source.lon_start) / source.lon_delta),
                        j_source_right = j_source_left + 1;
                    if (ambiguities[i * target.width + j] != nullptr) {
                        delete ambiguities[i * target.width + j];
                        ambiguities[i * target.width + j] = nullptr;
                        --number_of_amguities;
#ifdef DEBUG
                        fmt::print("Remove ambiguity records for {},{}\n", i, j);
#endif
                    };

                    T val[4] = { // lb lu rb ru
                        source.at(i_source_bottom, j_source_left),
                        source.at(i_source_upper, j_source_left),
                        source.at(i_source_bottom, j_source_right),
                        source.at(i_source_upper, j_source_right)
                    };
                    auto [lat, lon] = target.coordinate(i, j);
                    double coordinates[4][2];
                    std::tie(coordinates[0][0], coordinates[0][1]) = source.coordinate(i_source_bottom, j_source_left);
                    std::tie(coordinates[1][0], coordinates[1][1]) = source.coordinate(i_source_upper, j_source_left);
                    std::tie(coordinates[2][0], coordinates[2][1]) = source.coordinate(i_source_bottom, j_source_right);
                    std::tie(coordinates[3][0], coordinates[3][1]) = source.coordinate(i_source_upper, j_source_right);
                    switch (method)
                    {
                    case InterpolationMethod::Nearest:
                    {
                        double dist[4];
                        for (int k = 0; k < 4; ++k) {
                            dist[k] = (lat - coordinates[k][0]) * (lat - coordinates[k][0])
                                + (lon - coordinates[k][1]) * (lon - coordinates[k][1]);
                        }
                        size_t idx = 0;
                        for (int k = 1; k < 4; ++k) {
                            if (dist[k] < dist[idx]) {
                                idx = k;
                            }
                        }
                        target.data[target.width * i + j] = val[idx];
                        break;
                    }
                    // TODO: Add corresponding processing routine.
                    case InterpolationMethod::Bilinear:
                        break;
                    }

                    // Now we are quite certain about the value of this pixel.
                    mask[target.width * i + j] = 1;
                    --number_of_pixels_without_value;
                }
            };

            fmt::print("Processing outer extent.\n");
            size_t undertermined_points = 0;
            for (int i = bottom_out; i <= upper_out; ++i) {
                for (int j = left_out; j <= right_out; ++j) {
                    if (!(
                        ((i == bottom_out || i == upper_out) && (j >= left_out && j <= right_out)) ||
                        ((j == left_out || j == right_out) && (i >= bottom_out && i <= upper_out))
                        )) {
                        continue;
                    }

                    if (i < 0 || i >= target.height || j < 0 || j >= target.width) {
                        continue;
                    }

                    if (mask[target.width * i + j] == 1) {
                        continue;
                    };

                    auto [lat, lon] = target.coordinate(i, j);
                    int i_source_bottom = std::floor((target.lat_start + target.lat_delta * i - source.lat_start) / source.lat_delta),
                        i_source_upper = i_source_bottom + 1;
                    int j_source_left = std::floor((target.lon_start + target.lon_delta * j - source.lon_start) / source.lon_delta),
                        j_source_right = j_source_left + 1;
                    size_t valid_count = 0;
                    int coordinates[4][2] = {
                        i_source_bottom, j_source_left,
                        i_source_bottom, j_source_right,
                        i_source_upper, j_source_left,
                        i_source_upper, j_source_right
                    };
                    for (size_t k = 0; k < 4; ++k) {
                        if (coordinates[k][0] < 0 || coordinates[k][0] >= source.height
                            || coordinates[k][1] < 0 || coordinates[k][1] >= source.width) {
                            continue;
                        }
                        else {
                            Ambiguity* pa;
                            if (ambiguities[i * target.width + j] == nullptr) {
                                ++number_of_amguities;
                                pa = new Ambiguity;
                                ambiguities[i * target.width + j] = pa;
                            }

                            pa = ambiguities[i * target.width + j];
                            auto [s_lat, s_lon] = source.coordinate((size_t)coordinates[0], (size_t)coordinates[1]);
                            pa->candidates.push_back(std::make_tuple(s_lat, s_lon, source.at((size_t)coordinates[0], (size_t)coordinates[1])));
                            ++valid_count;
#ifdef DEBUG
                            fmt::print("For target point ({},{}), add point ({},{}) as a candidate.\n",
                                i, j, coordinates[k][0], coordinates[k][1]);
#endif
                        }
                    };
                    assert(valid_count <= 2);

                    undertermined_points++;
                };
            }
            fmt::print("Within this the extent of this .tif file, value for {} points in target raster are of ambiguity.\n", undertermined_points);
        };

        void resolve() {
            assert(number_of_amguities == 0);
            fmt::print("\n## Raster Merge Summary ##\nNumber of unresolved ambuities {}.\nNumber of pixels with no value {}.\n", number_of_amguities, number_of_pixels_without_value);
            return;
        };
    };
}