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

  auto doppler_getter = datatake.getNaiveDopplerGetter();
  double azimuth_interval =
      (datatake.azimuthEndTime - datatake.azimuthStartTime) /
      (datatake.numberOfAzimuthLines - 1);
  double range_interval =
      (datatake.rangeEndDistance - datatake.rangeStartDistance) /
      (datatake.numberOfRangePixels - 1);

  // 40: suggested in the article for ERS case, 26: range interval for ERS
  int range_polynomial_interval = 3 * 26 / range_interval,
      azimuth_polynomical_interval =
          5 * 30 / (azimuth_interval * datatake.states[0].vel.norm());

  range_polynomial_interval = 40;
  azimuth_polynomical_interval = 10;

  auto gen_grid_points = [](int end, int interval) {
    std::vector<int> grid_points;
    for (int i = 0; i < end; i += interval) {
      grid_points.push_back(i);
    }

    if (grid_points.back() != end - 1) {
      grid_points.push_back(end - 1);
    }

    return grid_points;
  };

  auto azimuth_grid_points = gen_grid_points(datatake.numberOfAzimuthLines,
                                             azimuth_polynomical_interval),
       range_grid_points = gen_grid_points(datatake.numberOfRangePixels,
                                           range_polynomial_interval);

  using Polynomial_type = Eigen::Matrix<double, 4, 2>;

  std::vector<Polynomial_type> polys(datatake.numberOfAzimuthLines *
                                     datatake.numberOfRangePixels);

  auto multi_polynomial =
      WHU::MultiPolynomial<WHU::StatePolynomial>::fromObjects(datatake.states,
                                                              6);

  double height_candidates[] = {1250, 3750, 6250, 8750};
  double height_offset = 5000, height_scale = 2000;

  Eigen::Matrix<double, 4, 4> tmp;
  for (int i = 0; i < 4; ++i) {
    double height_candidate = height_candidates[i];
    double effective_height = (height_candidate - height_offset) / height_scale,
           base = 1;
    for (int j = 3; j >= 0; --j) {
      tmp(i, j) = base;
      base *= effective_height;
    }
  }
  std::cout << tmp << std::endl;
  tmp = tmp.inverse().eval();

  // Gen polys
  for (auto azimuth_index : azimuth_grid_points) {
    fmt::print("Generating poly samples for line {}.\n", azimuth_index);
    auto azimuth_time =
        datatake.azimuthStartTime + datatake.azimuthLineTime() * azimuth_index;
    auto state = multi_polynomial.at(azimuth_time);

    for (auto range_index : range_grid_points) {
      double range = datatake.rangeStartDistance +
                     datatake.rangePixelSpacing() * range_index;

      Eigen::Matrix<double, 4, 2> Y;
      for (int i = 0; i < 4; ++i) {
        auto g = RDGeolocation(state, doppler_getter.doppler(range), range,
                               height_candidates[i]);
        Y(i, 0) = g.lon;
        Y(i, 1) = g.lat;
      }

      auto coeff = tmp * Y;

      polys[azimuth_index * datatake.numberOfRangePixels + range_index] = coeff;
    }
  };

  // Interpolating along lines
  for (auto azimuth_index : azimuth_grid_points) {
    fmt::print("Interpolating at azimuth line {}.\n", azimuth_index);
    for (int i = 0; i < range_grid_points.size() - 1; ++i) {
      int range_start = range_grid_points[i],
          range_end = range_grid_points[i + 1];
      auto left = polys[azimuth_index * datatake.numberOfRangePixels +
                        range_start],
           right =
               polys[azimuth_index * datatake.numberOfRangePixels + range_end];
      for (int j = range_start + 1; j < range_end; ++j) {
        double left_part = double(range_end - j) / (range_end - range_start),
               right_part = double(j - range_start) / (range_end - range_start);
        polys[azimuth_index * datatake.numberOfRangePixels + j] =
            left * left_part + right * right_part;
      }
    }
  };

  // Interpolating across lines
  auto interpolate_range_line = [&](int j) {
    if (j % 100 == 0) {
      fmt::print("Interpolating at range line {}.\n", j);
    }
    for (int i = 0; i < azimuth_grid_points.size() - 1; ++i) {
      int azimuth_start = azimuth_grid_points[i],
          azimuth_end = azimuth_grid_points[i + 1];

      auto up = polys[azimuth_start * datatake.numberOfRangePixels + j],
           down = polys[azimuth_end * datatake.numberOfRangePixels + j];
      for (int k = azimuth_start + 1; k < azimuth_end; ++k) {
        double up_part =
                   double(azimuth_end - k) / (azimuth_end - azimuth_start),
               down_part =
                   double(k - azimuth_start) / (azimuth_end - azimuth_start);
        polys[k * datatake.numberOfRangePixels + j] =
            up * up_part + down * down_part;
      }
    }
  };

  std::vector<int> indices;
  for (int i = 0; i < datatake.numberOfRangePixels; ++i) {
    indices.push_back(i);
  }
  std::for_each(std::execution::par_unseq, indices.begin(), indices.end(),
                [&](int i) { interpolate_range_line(i); });

  std::vector<double> height(datatake.numberOfAzimuthLines *
                             datatake.numberOfRangePixels, 0),
      azimuth_diff(datatake.numberOfAzimuthLines * datatake.numberOfRangePixels, 0),
      range_diff(datatake.numberOfAzimuthLines * datatake.numberOfRangePixels,
                   0);

  auto Eineder_geocode_line = [&](int i) {
    if (i % 100 == 0) {
      fmt::print("Eineder's geocoding for azimuth line {} finished.\n", i);
    }

    auto azimuth_time =
        datatake.azimuthStartTime + datatake.azimuthLineTime() * i;
    auto state = multi_polynomial.at(azimuth_time);

    for (int j = 0; j < datatake.numberOfRangePixels; ++j) {
      auto poly = polys[i * datatake.numberOfRangePixels + j];
      double range =
          datatake.rangeStartDistance + datatake.rangePixelSpacing() * j;
      // fmt::print("Poly at {},{}: \n", i, j);
      // std::cout << poly << std::endl;

      constexpr double convergent_condition = 1e-5;
      Eigen::Matrix<double, 1, 4> coeffs;

      Eigen::Matrix<double, 1, 2> horizontal_coordinates;
      double height_from_dem;

      // Bi-section, monotonically decreasing
      auto get_diff = [&](double height) {
        double base = 1;
        double height_tmp = (height - height_offset) / height_scale;
        for (int k = 3; k >= 0; --k) {
          coeffs(0, k) = base;
          base *= height_tmp;
        };

        // lon lat
        horizontal_coordinates = (coeffs * poly).eval();

        height_from_dem = dem.at_bspline(horizontal_coordinates(0, 1),
                                            horizontal_coordinates(0, 0));

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

      GeodeticCoordinate geo;
      geo.lon = horizontal_coordinates(0, 0);
      geo.lat = horizontal_coordinates(0, 1);
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

  indices.resize(0);
  for (int i = 0; i < datatake.numberOfAzimuthLines; ++i) {
    indices.push_back(i);
  }
  std::for_each(std::execution::par_unseq, indices.begin(), indices.end(),
                [&](int i) { Eineder_geocode_line(i); });

  timer.stop();
  fmt::print("Geocoding {} ms.\n", timer.elapsedMilliseconds());

  std::ofstream height_file((output_dir / "eineder_height.bin").c_str(),
                            std::ios::binary | std::ios::out);
  height_file.write((char*)(height.data()), sizeof(double) * height.size());
  height_file.close();

  std::ofstream azimuth_diff_file((output_dir / "eineder_azimuth_diff.bin").c_str(),
                            std::ios::binary | std::ios::out);
  azimuth_diff_file.write((char*)(azimuth_diff.data()),
                    sizeof(double) * azimuth_diff.size());
  azimuth_diff_file.close();

  std::ofstream range_diff_file((output_dir / "eineder_range_diff.bin").c_str(),
                            std::ios::binary | std::ios::out);
  range_diff_file.write((char*)(range_diff.data()),
                    sizeof(double) * range_diff.size());
  range_diff_file.close();


  fmt::print("Azimutn lines: {}, range samples: {}.\n",
             datatake.numberOfAzimuthLines, datatake.numberOfRangePixels);
  return 0;
};