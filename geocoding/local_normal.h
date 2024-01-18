#pragma once

#include <array>
#include <vector>

#include "th_datatake.h"
#include "th_dem.h"
#include "th_interp.h"

struct LocalNormal {
  LocalNormal() = default;

  LocalNormal split(const LocalNormal& other, double ratio) {
    LocalNormal tmp;
    for (int i = 0; i < 8; ++i) {
      tmp.data[i] = data[i] + (other.data[i] - data[i]) * ratio;
    }
    return tmp;
  }

  LocalNormal operator+(const LocalNormal& other) {
    LocalNormal tmp;
    for (int i = 0; i < 8; ++i) {
      tmp.data[i] = data[i] + other.data[i];
    }
    return tmp;
  }

  LocalNormal operator*(const LocalNormal& other) {
    LocalNormal tmp;
    for (int i = 0; i < 8; ++i) {
      tmp.data[i] = data[i] * other.data[i];
    }
    return tmp;
  }

  LocalNormal operator/(const LocalNormal& other) {
    LocalNormal tmp;
    for (int i = 0; i < 8; ++i) {
      tmp.data[i] = data[i] / other.data[i];
    }
    return tmp;
  }

  LocalNormal operator-(const LocalNormal& other) {
    LocalNormal tmp;
    for (int i = 0; i < 8; ++i) {
      tmp.data[i] = data[i] - other.data[i];
    }
    return tmp;
  }

  template <size_t N>
  static LocalNormal from_array(const std::array<double, N>& arr) {
    static_assert(N >= 6);

    LocalNormal tmp;
    for (size_t i = 0; i < std::min<size_t>(N, 8); ++i) {
      tmp.data[i] = arr[i];
    }

    return tmp;
  }

  // as a hint for auto vectorizer,
  // padding to 512 bits to avoid overlapping between adjacent elements in
  // array.
  // __align might work
  double data[8];

  double hori_lon() const { return data[0]; };

  double hori_lat() const { return data[1]; };

  double hgt_lon() const { return data[2]; };

  double hgt_lat() const { return data[3]; };

  double orbit_lon() const { return data[4]; };

  double orbit_lat() const { return data[5]; };
};

// More sophisticated interpolation method to be involved?
// Retrive local normal approximation by (azimuth_line_number, slant_range)
// (size_t, double)
struct LocalNormalGrid {
  size_t num_range_samples;
  double range_start, range_end, range_step;

  std::vector<LocalNormal> normals;

  LocalNormal at(size_t line, double range) {
    int left = (range - range_start) / range_step, right = left + 1;
    double margin = (range - range_start) / range_step - left;

    if (left < 0) {
      return normals[line * num_range_samples];
    }
    if (right >= num_range_samples) {
      return normals[(line + 1) * num_range_samples - 1];
    }

    return normals[line * num_range_samples + left].split(
        normals[line * num_range_samples + right], margin);
  };

  LocalNormalGrid(const WHU::Datatake& datatake, double range_margin,
                  size_t num_azimuth_effective_samples,
                  size_t num_range_samples, WHU::DEM& dem)
      : num_range_samples(num_range_samples),
        range_start(datatake.rangeStartDistance - range_margin),
        range_end(datatake.rangeEndDistance + range_margin),
        range_step((range_end - range_start) / (num_range_samples - 1)),
        normals(datatake.numberOfAzimuthLines * num_range_samples) {
    datatake.states[0].time.assert_same_day(datatake.states.back().time);

    auto multi_polynomial =
        WHU::MultiPolynomial<WHU::StatePolynomial>::fromObjects(datatake.states,
                                                                6);

    double slant_range_delta = 1e-2, height_delta = 1e-2, orbit_delta = 1e-7;

    struct TimeNormalTie {
      double time;
      LocalNormal normal;

      static auto getTimeAccessor() {
        auto tmp = [](const TimeNormalTie& instance) { return instance.time; };
        return tmp;
      };

      static auto getDataAccessors() {
        auto get_1 = [](const TimeNormalTie& instance) {
          return instance.normal.data[0];
        };
        auto get_2 = [](const TimeNormalTie& instance) {
          return instance.normal.data[1];
        };
        auto get_3 = [](const TimeNormalTie& instance) {
          return instance.normal.data[2];
        };
        auto get_4 = [](const TimeNormalTie& instance) {
          return instance.normal.data[3];
        };
        auto get_5 = [](const TimeNormalTie& instance) {
          return instance.normal.data[4];
        };
        auto get_6 = [](const TimeNormalTie& instance) {
          return instance.normal.data[5];
        };
        std::array<std::function<double(const TimeNormalTie&)>, 6> tmp = {
            get_1, get_2, get_3, get_4, get_5, get_6};
        return tmp;
      };
    };

    // Sparse normal samples in azimuth direction.
    std::vector<std::vector<TimeNormalTie>> normal_samples(
        num_range_samples,
        std::vector<TimeNormalTie>(num_azimuth_effective_samples));

    auto doppler_getter = datatake.getNaiveDopplerGetter();

    auto gen_sample_line = [&](size_t i) {
      double azimuth_delta =
          (datatake.azimuthEndTime - datatake.azimuthStartTime) /
          (num_azimuth_effective_samples - 1);
      auto azimuth_curr = datatake.azimuthStartTime + azimuth_delta * i;

      auto state = multi_polynomial.at(azimuth_curr);
      auto state_delta = multi_polynomial.at(state.time + orbit_delta);
      double range_curr = range_start,
             range_delta = (range_end - range_start) / (num_range_samples - 1);

      auto get_geodetic_pos = [&](double rho) {
        // A very loose one.
        constexpr double convergent_condition = 1;

        double height_from_dem;
        WHU::GeodeticCoordinate geo_iter;

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

      for (int j = 0; j < num_range_samples; ++j) {
        auto g1 = get_geodetic_pos(range_curr);
        if (g1.height == 0) {
          std::cout << "Not converged at " << i << '\t' << j << std::endl;
          LocalNormal tmp;
          tmp.data[0] = std::nan("");

          TimeNormalTie tie;
          tie.time = azimuth_curr.seconds;
          tie.normal = tmp;

          normal_samples[j][i] = tie;

          range_curr += range_delta;
          continue;
        }

        // p1: bottom_left p2: bottom_right p3: top
        auto p1 = RDGeolocation(
                 state,
                 doppler_getter.doppler(range_curr - slant_range_delta / 2),
                 range_curr - slant_range_delta / 2,
                 g1.height - height_delta / 2),
             p2 = RDGeolocation(
                 state,
                 doppler_getter.doppler(range_curr + slant_range_delta / 2),
                 range_curr + slant_range_delta / 2,
                 g1.height - height_delta / 2),
             p3 = RDGeolocation(state, doppler_getter.doppler(range_curr),
                                range_curr, g1.height + height_delta / 2),
             p4 = RDGeolocation(state, doppler_getter.doppler(range_curr),
                                range_curr, g1.height),
             p5 = RDGeolocation(
                 state_delta,
                 doppler_getter.doppler(range_curr),
                 range_curr, g1.height);

        auto hori_delta_lon = p2.lon - p1.lon, hori_delta_lat = p2.lat - p1.lat,
             slant_delta_lon = p3.lon - p1.lon,
             slant_delta_lat = p3.lat - p1.lat;
        double k =
            (hori_delta_lon * slant_delta_lon +
             hori_delta_lat * slant_delta_lat) /
            (hori_delta_lon * hori_delta_lon + hori_delta_lat * hori_delta_lat);
        double verti_delta_lon = slant_delta_lon - k * hori_delta_lon,
               verti_delta_lat = slant_delta_lat - k * hori_delta_lat;
        double test =
            hori_delta_lon * verti_delta_lon + hori_delta_lat * verti_delta_lat;

        double orbit_delta_lon = p5.lon - p4.lon,
               orbit_delta_lat = p5.lat - p4.lat;

        LocalNormal tmp;
        tmp.data[0] = hori_delta_lon / slant_range_delta;
        tmp.data[1] = hori_delta_lat / slant_range_delta;
        tmp.data[2] = verti_delta_lon / height_delta;
        tmp.data[3] = verti_delta_lat / height_delta;
        tmp.data[4] = orbit_delta_lon / orbit_delta;
        tmp.data[5] = orbit_delta_lat / orbit_delta;

#define N_6CURIOUS
#ifdef CURIOUS
        for (int j = 0; j < 20; ++j) {
          auto p1 = RDGeolocation(state, doppler_getter.doppler(range_curr),
                                  range_curr, g1.height + (j - 10) * 100);
          for (int i = 0; i < 20; ++i) {
            double height_delta = i * 100;
            auto test_tmp = p1;
            test_tmp.lon += tmp.data[2] * height_delta;
            test_tmp.lat += tmp.data[3] * height_delta;
            test_tmp.height += height_delta;
            // auto t = state.estimateTimeDelta(geodetic2cartesian(test_tmp));
            auto pos_tmp = geodetic2cartesian(test_tmp);
            auto t = state.estimateTimeDeltaWithDoppler(
                pos_tmp, doppler_getter.doppler(pos_tmp.distance2(state.pos)));
            fmt::print(
                "Origin deviation : {}, Height delta : {}, azimuth_delta {}.\n",
                (j - 10) * 100, height_delta, t * 1e4);
          }
        }
#endif

        TimeNormalTie tie;
        tie.time = azimuth_curr.seconds;
        tie.normal = tmp;

        normal_samples[j][i] = tie;

        range_curr += range_delta;
      };

      azimuth_curr += azimuth_delta;
    };

    std::vector<size_t> indices(num_azimuth_effective_samples);
    for (size_t i = 0; i < num_azimuth_effective_samples; ++i) {
      indices[i] = i;
    }
    std::for_each(std::execution::par_unseq, indices.begin(), indices.end(),
                  [&](size_t i) { gen_sample_line(i); });

    static constexpr size_t interp_degree = 3;

    auto interp_range_line = [&](size_t j) {
      auto multipol = WHU::MultiPolynomial<WHU::PolynomialNT<6>>::fromObjects(
          normal_samples[j], interp_degree + 1);
      double azimuth_curr = datatake.azimuthStartTime.seconds,
             azimuth_delta =
                 (datatake.azimuthEndTime - datatake.azimuthStartTime) /
                 (datatake.numberOfAzimuthLines - 1);

      for (size_t i = 0; i < datatake.numberOfAzimuthLines; ++i) {
        auto tmp = LocalNormal::from_array(multipol.at(azimuth_curr));
        azimuth_curr += azimuth_delta;

        normals[i * num_range_samples + j] = tmp;
      };
    };

    indices.resize(num_range_samples);
    for (size_t i = 0; i < num_range_samples; ++i) {
      indices[i] = i;
    }

    std::for_each(std::execution::par_unseq, indices.begin(), indices.end(),
                  [&](size_t j) { interp_range_line(j); });
  };
};
