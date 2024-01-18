#include <fmt/format.h>
#include <th_constants.h>
#include <th_datatake.h>
#include <th_dem.h>
#include <th_geometric.h>
#include <th_interp.h>
#include <th_raster.h>
#include <th_statepolynomial.h>
#include <th_timer.h>

#include <Eigen/Dense>
#include <array>
#include <cxxopts.hpp>
#include <execution>
#include <fstream>
#include <memory>
#include <thread>

#include "local_normal.h"

// This file is too ugly!

using namespace WHU;

#ifdef DEBUG
#undef DEBUG
#endif
#ifdef NDEBUG
#undef NDEBUG
#endif

#define NDEBUG

// Except for the output format, nothing is modified.
int main(int argc, char** argv) {
  cxxopts::Options options("Geocoding",
                           "A fast and accurate forward geocoding method");
  options.add_options()(
      "d,datatake",
      "datatake file containing neccessary geometric information.",
      cxxopts::value<fs::path>())(
      "e,dem", "DEM covering the target ground scene.",
      cxxopts::value<
          fs::path>())("o,output", "Directory containing output binaries.",
                       cxxopts::value<
                           fs::path>())("r,override",
                                        "Enable overriding existing "
                                        "files.")("s, single",
                                                  "Disable multithreading.",
                                                  cxxopts::value<bool>());

  if (argc == 1) {
    std::cout << options.help();
    return 0;
  }
  auto arguments = options.parse(argc, argv);

  auto input_datatake = arguments["d"].as<fs::path>();
  auto input_dem = arguments["e"].as<fs::path>();
  auto output_dir = arguments["o"].as<fs::path>();
  auto enabling_override = arguments["r"].as<bool>();
  auto disable_multithreading = arguments["s"].as<bool>();

  if (!fs::is_regular_file(input_datatake)) {
    fmt::print("Cannot find input datatake file {}.\n",
               input_datatake.string());
    return 1;
  }

  if (!fs::is_regular_file(input_dem)) {
    fmt::print("Cannot find input dem file {}.\n", input_dem.string());
    return 1;
  }

  if (fs::is_regular_file(output_dir) ||
      (fs::is_directory(output_dir) && !fs::is_empty(output_dir) &&
       !enabling_override)) {
    fmt::print("Output dir {} is a file or a directory that is not empty.\n",
               input_dem.string());
    return 1;
  }

  if (!fs::is_directory(output_dir)) {
    fs::create_directory(output_dir);
  };

  Timer timer;

  bool longitude_from_zero = false;

  // Due to unrecognized problems in GAMMA, zero values should be
  static constexpr double zero_height_modification = 0.0001;

  timer.start();
  auto dem = DEM::fromBoostBinary(input_dem.string());
  for (auto& h : dem.data) {
    if (h == 0.) {
      h = zero_height_modification;
    }
  };
  // dem.set_not_use_spline();
  dem.init_bspline();

  if (dem.longitudeEndArc > 180 * 3600) {
    longitude_from_zero = true;
    set_longitude_from_zero(true);
  };
  timer.stop();
  fmt::print("Importing DEM costs {} ms.\n", timer.elapsedMilliseconds());

  auto datatake = Datatake::fromBoostText(input_datatake.string());
  // datatake.numberOfRangePixels *= 3;
  // datatake.numberOfAzimuthLines *= 4;

  // datatake.numberOfAzimuthLines /= 40;

  // Debug
  // datatake.numberOfAzimuthLines = 30;

  auto doppler_getter = datatake.getNaiveDopplerGetter();
  double azimuth_interval =
      (datatake.azimuthEndTime - datatake.azimuthStartTime) /
      (datatake.numberOfAzimuthLines - 1);
  double range_interval =
      (datatake.rangeEndDistance - datatake.rangeStartDistance) /
      (datatake.numberOfRangePixels - 1);

  std::vector<double> height(datatake.numberOfAzimuthLines *
                             datatake.numberOfRangePixels, 0),
      azimuth_diff(datatake.numberOfAzimuthLines * datatake.numberOfRangePixels, 0),
      range_diff(datatake.numberOfAzimuthLines * datatake.numberOfRangePixels,
                   0);

  auto multi_polynomial =
      WHU::MultiPolynomial<WHU::StatePolynomial>::fromObjects(datatake.states,
                                                              6);

  auto Eineder_geocode_line = [&](int i) {
    if (i % 100 == 0) {
      fmt::print("Eineder's geocoding for azimuth line {} finished.\n", i);
    }

    auto azimuth_time =
        datatake.azimuthStartTime + datatake.azimuthLineTime() * i;
    auto state = multi_polynomial.at(azimuth_time);

    for (int j = 0; j < datatake.numberOfRangePixels; ++j) {
      double range =
          datatake.rangeStartDistance + datatake.rangePixelSpacing() * j;
      // fmt::print("Poly at {},{}: \n", i, j);
      // std::cout << poly << std::endl;

      constexpr double convergent_condition = 1e-3;

      double height_from_dem;
      GeodeticCoordinate geo_iter;

      // Bi-section, monotonically decreasing
      auto get_diff = [&](double height) {
        geo_iter = RDGeolocation(state, 0, range, height, convergent_condition);
        
        height_from_dem = dem.at_bspline(geo_iter.lat, geo_iter.lon);

        return height_from_dem - height;
      };

      double left = -500, right = 10000;
      while (right - left > convergent_condition) {
        double mid = (right + left) / 2;
        if (get_diff(mid) < 0) {
          right = mid;
        } else {
          left = mid;
        }
      }

      /*
      bool found = true;
      while (std::abs(height_now - height_prev) > convergent_condition) {
        height_prev = height_now;


        if (++iter_count == 200) {
          found = false;
          break;
        }
      };*/

      GeodeticCoordinate geo;
      geo.lon = geo_iter.lon;
      geo.lat = geo_iter.lat;
      geo.height = height_from_dem;

      auto Cartesian = geodetic2cartesian(geo);

      double sample_azimuth_diff =
          state.estimateTimeDelta(Cartesian) * state.vel.norm();
      double sample_range_diff = Cartesian.distance2(state.pos) - range;

      /*
      fmt::print("Azimuth diff: {}, range diff: {}, height: {}.\n",
                 sample_azimuth_diff, sample_range_diff, height_from_dem);
      fmt::print("Converged at {},{}.\n", i, j);*/

      height[i * datatake.numberOfRangePixels + j] = height_from_dem;
      azimuth_diff[i * datatake.numberOfRangePixels + j] = sample_azimuth_diff;
      range_diff[i * datatake.numberOfRangePixels + j] = sample_range_diff;
    }
  };

  timer.start();

  std::vector<int> indices;
  indices.resize(0);
  for (int i = 0; i < datatake.numberOfAzimuthLines; ++i) {
    indices.push_back(i);
  }
  std::for_each(std::execution::par_unseq, indices.begin(), indices.end(),
                [&](int i) { Eineder_geocode_line(i); });

  timer.stop();
  fmt::print("Geocoding {} ms.\n", timer.elapsedMilliseconds());

  std::ofstream height_file((output_dir / "plain_height.bin").c_str(),
                            std::ios::binary | std::ios::out);
  height_file.write((char*)(height.data()), sizeof(double) * height.size());
  height_file.close();

  std::ofstream azimuth_diff_file((output_dir / "plain_azimuth_diff.bin").c_str(),
                            std::ios::binary | std::ios::out);
  azimuth_diff_file.write((char*)(azimuth_diff.data()),
                    sizeof(double) * azimuth_diff.size());
  azimuth_diff_file.close();

  std::ofstream range_diff_file((output_dir / "plain_range_diff.bin").c_str(),
                            std::ios::binary | std::ios::out);
  range_diff_file.write((char*)(range_diff.data()),
                    sizeof(double) * range_diff.size());
  range_diff_file.close();

  fmt::print("Azimutn lines: {}, range samples: {}.\n",
             datatake.numberOfAzimuthLines, datatake.numberOfRangePixels);
  return 0;
};