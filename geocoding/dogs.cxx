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

  timer.start();
  LocalNormalGrid grid(datatake, 1e3, 40, 40, dem);
  timer.stop();
  fmt::print("Generating normal grid costs {} ms. \n",
             timer.elapsedMilliseconds());

  auto t = MultiPolynomial<StatePolynomial>::fromObjects(datatake.states, 6);

  // Baseline vector for the detection of shadow region.
  auto azimuth_mid_time =
      datatake.azimuthStartTime +
      (datatake.azimuthEndTime - datatake.azimuthStartTime) / 2;
  auto range_mid_dist =
      (datatake.rangeStartDistance + datatake.rangeEndDistance) / 2;
  auto midpoint_state = t.at(azimuth_mid_time);
  auto base_vec = geodetic2cartesian(
      RDGeolocation(midpoint_state, doppler_getter.doppler(range_mid_dist),
                    range_mid_dist, 800000));
  base_vec = base_vec - midpoint_state.pos;
  base_vec = base_vec * 1 / base_vec.norm();

  std::vector<GeodeticCoordinate> geocode_results(
      datatake.numberOfAzimuthLines * datatake.numberOfRangePixels);

  // Call to memset may be optimized out!
  std::memset(geocode_results.data(), 0,
              sizeof(GeodeticCoordinate) * geocode_results.size());

  struct DuplicatedPoint {
    int azimuth;
    int range;
    GeodeticCoordinate geo;
  };

  std::vector<std::vector<DuplicatedPoint>> duplicated(
      datatake.numberOfAzimuthLines);
  std::vector<std::vector<DuplicatedPoint>> out_of_boundary(
      datatake.numberOfAzimuthLines);
  std::vector<std::vector<std::array<double, 2>>> layovers(
      datatake.numberOfAzimuthLines);

  
  std::vector<double> height(
      datatake.numberOfAzimuthLines * datatake.numberOfRangePixels, 0),
      azimuth_diff(datatake.numberOfAzimuthLines * datatake.numberOfRangePixels,
                   0),
      range_diff(datatake.numberOfAzimuthLines * datatake.numberOfRangePixels,
                 0);

  // Geocode one line very fast, line: line number; ratio: upscaling ratio;
  // init_status_points: number of points used to init estimator
  auto geocode_one_line = [&](size_t line, double range_base_interval,
                              double deviation_maximum = 0.01) {
    auto state = t.at(datatake.azimuthStartTime + azimuth_interval * line);

    // GeodeticCoordinate tmp_debug{0, 0, 0};

    auto get_geodetic_pos = [&](double rho) {
      constexpr double convergent_condition = 1e-5;

      double height_from_dem;
      GeodeticCoordinate geo_iter;

      // Bi-section, monotonically decreasing
      auto get_diff = [&](double height) {
        geo_iter = RDGeolocation(state, 0, rho, height);

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

      geo_iter.height = height_from_dem;
      return geo_iter;
    };

    // Some data pixels may appear multiple times on the ground.
    std::vector<GeodeticCoordinate> slant_range_points(
        datatake.numberOfRangePixels);
    std::vector<bool> geocoded(datatake.numberOfRangePixels, false);
    int num_geocoded = 0, num_geocoded_repeated = 0;

    // Hyper parameters
    int margin = 500, double_step_after_consecutive_valid_samples = 5;
    double ignore_insignificant_layover = 0.45;
    double step_ratio_lower_bound = 1;
    while (step_ratio_lower_bound * range_base_interval >
           ignore_insignificant_layover) {
      step_ratio_lower_bound /= 2;
    }

    bool found = false;
    double curr_range, curr_height, curr_index;
    GeodeticCoordinate curr_geo;

    // Never mind, when stupid search method was adopted.
    int find_init_sample_record = 0;
    while (!found) {
      if (++find_init_sample_record % 5 == 0) {
        fmt::print("{}-th time to get the first sample of line {}.\n",
                   find_init_sample_record, line);
      }
      curr_range = datatake.rangeStartDistance - range_interval * margin;
      curr_geo = get_geodetic_pos(curr_range);
      curr_index = -margin;
      if (curr_geo.height) {
        curr_height = dem.at_bspline(curr_geo.lat, curr_geo.lon);
        found = true;
      } else {
        // fmt::print("Retry after all finding algorithms failed for line
        // {}.\n", line);
        margin -= 10;
      }
    }

    // This term should be always decreasing.
    double look_angle_cosine;
    bool not_first_sample_out_of_shadow = false;
    {
      auto tmp = geodetic2cartesian(curr_geo);
      auto tmp_vec = tmp - state.pos;
      look_angle_cosine = base_vec.dot(tmp_vec) / tmp_vec.norm();
    }

    // The slope
    double curr_range_deviation_ratio =
               0,  // expected_range_deviation_per_unit_range_delta
        step_ratio = 1;

    // Penalty and bonus for ground sampling interval. 
    int consecutive_valid_count = 0;
    auto step_ratio_penalty = [&]() {
      double tmp = step_ratio / 2;
      if (tmp < step_ratio_lower_bound) {
        return false;
      } else {
        step_ratio = tmp;
        return true;
      }
    };
    auto step_ratio_boost = [&](bool good) {
      if (good) {
        ++consecutive_valid_count;
      } else {
        consecutive_valid_count = 0;
        return;
      }

      if (consecutive_valid_count >=
              double_step_after_consecutive_valid_samples &&
          step_ratio <= 0.5) {
        consecutive_valid_count = 0;
        step_ratio *= 2;
      };
    };

    bool in_layover = false;
    double into_layover_range, minimum_layover_range;

    int next_range_pixel_index = 0;

    assert(curr_range < datatake.rangeStartDistance);

    int iter_count = 0;
    // Main loop
    while (curr_range < datatake.rangeEndDistance) {
      ++iter_count;
      auto normal = grid.at(line, curr_range);

      double expected_range_delta = range_base_interval * step_ratio;
      double expected_range_deviation =
          expected_range_delta * curr_range_deviation_ratio;

      // Move along the horizontal normal
      auto next_geo = curr_geo;
      next_geo.lon += normal.hori_lon() * expected_range_delta;
      next_geo.lat += normal.hori_lat() * expected_range_delta;

      // Move along the height-variation-induced target-moving normal.
      double height_delta =
          dem.at_bspline(next_geo.lat, next_geo.lon) - curr_height;
      next_geo.lon += height_delta * normal.hgt_lon();
      next_geo.lat += height_delta * normal.hgt_lat();

      // Move along the normal perpendicular to local Dopper plane.
      next_geo.height = dem.at_bspline(next_geo.lat, next_geo.lon);
      auto pos = geodetic2cartesian(next_geo);
      auto residual = state.estimateTimeDeltaWithDoppler(
          pos, doppler_getter.doppler(pos.distance2(state.pos)));
      next_geo.lon -= residual * normal.orbit_lon();
      next_geo.lat -= residual * normal.orbit_lat();
      next_geo.height = dem.at_bspline(next_geo.lat, next_geo.lon);
      auto next_pos = geodetic2cartesian(next_geo);

      auto next_look_vector = next_pos - state.pos;
      double next_range = next_look_vector.norm(),
             delta_range = next_range - curr_range;
      double next_range_deviation_ratio =
          (expected_range_delta - delta_range) / expected_range_delta;

      if (std::abs(next_range_deviation_ratio - curr_range_deviation_ratio) >
              deviation_maximum &&
          step_ratio_penalty()) {
        continue;
      }
      // Now we are confident that this sample is valid.

      step_ratio_boost(
          std::abs(next_range_deviation_ratio - curr_range_deviation_ratio) <
          deviation_maximum / 2);

      // But it might be in a shadow.
      double next_index =
          (next_range - datatake.rangeStartDistance) / range_interval;
      auto tmp_look_angle_cosine = next_look_vector.dot(base_vec) / next_range;

      if (tmp_look_angle_cosine > look_angle_cosine) {
        look_angle_cosine = tmp_look_angle_cosine;

        if (in_layover) {
          minimum_layover_range = std::min(minimum_layover_range, next_range);
        }

        // Just move into a layover region.
        if (delta_range < 0 && !in_layover) {
          in_layover = true;
          into_layover_range = curr_range;
          minimum_layover_range = curr_range;
        }

        // Just move out of a layover region.
        if (delta_range > 0 && in_layover) {
          in_layover = false;
          layovers[line].push_back({minimum_layover_range, into_layover_range});
        };

        if (not_first_sample_out_of_shadow) {
          // Interpolate to get samples in Slant-range coordinate system.
          // next_index = (next_range - datatake.rangeStartDistance) /
          // range_interval;
          double index_lower_d = curr_index, index_upper_d = next_index;
          bool swapped = 0;
          if (index_lower_d > index_upper_d) {
            std::swap(index_lower_d, index_upper_d);
            swapped = 1;
          }
          int index_lower = std::ceil(index_lower_d),
              index_upper = std::floor(index_upper_d);
          for (int i = index_lower; i <= index_upper; ++i) {
            double slant_range =
                datatake.rangeStartDistance + i * range_interval;
            double curr_coeff =
                       std::abs((slant_range - next_range) / delta_range),
                   next_coeff =
                       std::abs((slant_range - curr_range) / delta_range);
            GeodeticCoordinate tmp;
            tmp.lat = curr_coeff * curr_geo.lat + next_coeff * next_geo.lat;
            tmp.lon = curr_coeff * curr_geo.lon + next_coeff * next_geo.lon;
            tmp.height = dem.at_bspline(tmp.lat, tmp.lon);

            if (i < 0 || i >= datatake.numberOfRangePixels) {
              DuplicatedPoint point{line, -1, tmp};
              out_of_boundary[line].push_back(point);
              continue;
            }

            if (geocoded[i]) {
              ++num_geocoded_repeated;
              DuplicatedPoint point{line, i, tmp};
              duplicated[line].push_back(point);
            } else {
              geocode_results[line * datatake.numberOfRangePixels + i] = tmp;

              height[line * datatake.numberOfRangePixels + i] = tmp.height;
              auto tmpxyz = geodetic2cartesian(tmp);
              double diff_time = state.estimateTimeDeltaWithDoppler(
                         tmpxyz,
                         doppler_getter.doppler(tmpxyz.distance2(state.pos))) * state.vel.norm(),
                     diff_range =
                         tmpxyz.distance2(state.pos) -
                         (datatake.rangeStartDistance + i * range_interval);
              azimuth_diff[line * datatake.numberOfRangePixels + i] = diff_time;
              range_diff[line * datatake.numberOfRangePixels + i] = diff_range;

              geocoded[i] = true;
              ++num_geocoded;
            }
          }
        }
        not_first_sample_out_of_shadow = true;
      } else {
        not_first_sample_out_of_shadow = false;
      }

      // Update parameters.
      curr_range = next_range;
      curr_geo = next_geo;
      curr_height = next_geo.height;
      curr_range_deviation_ratio = next_range_deviation_ratio;
      curr_index = next_index;
      // samples.push_back(curr_geo);
    };
  };

  timer.start();
  std::vector<size_t> indices;
  std::atomic_size_t line_counter = 0;
  for (size_t i = 0; i < datatake.numberOfAzimuthLines; ++i) {
    indices.push_back(i);
  }

  if (disable_multithreading) {
    std::for_each(
        std::execution::seq, indices.begin(), indices.end(), [&](size_t index) {
          double base_spacing = datatake.rangePixelSpacing();
          if (base_spacing < 3) {
            base_spacing = 3;
          }
          geocode_one_line(index, base_spacing, 0.01);
          // dummy += geocode_one_line(index, 0.3);

          line_counter++;
          if (line_counter % 100 == 0) {
            fmt::print("{} / {} lines geocoded.\n", double(line_counter),
                       datatake.numberOfAzimuthLines);
          }
        });
  } else {
    std::for_each(std::execution::par_unseq, indices.begin(), indices.end(),
                  [&](size_t index) {
                    double base_spacing = datatake.rangePixelSpacing();
                    if (base_spacing < 3) {
                      base_spacing = 3;
                    }
                    geocode_one_line(index, base_spacing, 0.01);

                    line_counter++;
                    if (line_counter % 100 == 0) {
                      fmt::print("{} / {} lines geocoded.\n",
                                 double(line_counter),
                                 datatake.numberOfAzimuthLines);
                    }
                  });
  }
  timer.stop();

  double samples_per_line =
      geocode_results.size() / datatake.numberOfAzimuthLines;
  fmt::print("Get {} samples for {} lines used {} milliseconds. Average : {}\n",
             geocode_results.size(), datatake.numberOfAzimuthLines,
             timer.elapsedMilliseconds(), samples_per_line);

  std::vector<unsigned char> is_multiple_points(
      datatake.numberOfAzimuthLines * datatake.numberOfRangePixels, false);
  for (auto& tmp : duplicated) {
    for (auto& p : tmp) {
      is_multiple_points[p.azimuth * datatake.numberOfRangePixels + p.range] =
          true;
    }
  }

  std::ofstream multiple_file((output_dir / "dog_multiple.bin").c_str(),
                            std::ios::binary | std::ios::out);
  multiple_file.write((char*)(is_multiple_points.data()),
                    sizeof(unsigned char) * height.size());
  multiple_file.close();


  std::ofstream height_file((output_dir / "dog_height.bin").c_str(),
                            std::ios::binary | std::ios::out);
  height_file.write((char*)(height.data()), sizeof(double) * height.size());
  height_file.close();

  std::ofstream azimuth_diff_file(
      (output_dir / "dog_azimuth_diff.bin").c_str(),
      std::ios::binary | std::ios::out);
  azimuth_diff_file.write((char*)(azimuth_diff.data()),
                          sizeof(double) * azimuth_diff.size());
  azimuth_diff_file.close();

  std::ofstream range_diff_file((output_dir / "dog_range_diff.bin").c_str(),
                                std::ios::binary | std::ios::out);
  range_diff_file.write((char*)(range_diff.data()),
                        sizeof(double) * range_diff.size());
  range_diff_file.close();

  fmt::print("Azimutn lines: {}, range samples: {}.\n",
             datatake.numberOfAzimuthLines, datatake.numberOfRangePixels);


  return 0;
};