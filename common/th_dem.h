#pragma once

#include <cmath>
#include <optional>

#include <fmt/format.h>
#include <Eigen/Dense>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include "th_bspline.h"
#include "th_serialize.h"

namespace WHU {
struct DEM : public ISerializable<DEM> {
  int longitudeStartArc;
  int longitudeEndArc;
  int latitudeStartArc;
  int latitudeEndArc;

  size_t numberLongitudePixels;
  size_t numberLatitudePixels;

  bool not_use_spline = false;

  std::vector<float> data;

  std::optional<BSpline2D> spline_interp;

  void put(int lat_arc, int lon_arc, float val) {
    if (lat_arc < latitudeStartArc || lat_arc > latitudeEndArc ||
        lon_arc < longitudeStartArc || lon_arc > longitudeEndArc) {
      return;
    };

    data[(lat_arc - latitudeStartArc) * numberLongitudePixels +
         (lon_arc - longitudeStartArc)] = val;
  }

  double at(double lat, double lon) const {
    lat *= 3600, lon *= 3600;
    if (lat < latitudeStartArc || lat > latitudeEndArc ||
        lon < longitudeStartArc || lon > longitudeEndArc) {
      return std::nan("");
    }

    int lat_lo = std::floor(lat), lon_lo = std::floor(lon);
    double lat_dl = lat - lat_lo, lat_dh = 1 - lat_dl, lon_dl = lon - lon_lo,
           lon_dh = 1 - lon_dl;
    double lb = data[(lat_lo - latitudeStartArc) * numberLongitudePixels +
                     (lon_lo - longitudeStartArc)],
           lh = data[(lat_lo + 1 - latitudeStartArc) * numberLongitudePixels +
                     (lon_lo - longitudeStartArc)],
           rb = data[(lat_lo - latitudeStartArc) * numberLongitudePixels +
                     (lon_lo + 1 - longitudeStartArc)],
           rh = data[(lat_lo + 1 - latitudeStartArc) * numberLongitudePixels +
                     (lon_lo + 1 - longitudeStartArc)];
    return lb * lat_dh * lon_dh + lh * lat_dl * lon_dh + rb * lat_dh * lon_dl +
           rh * lat_dl * lon_dl;
  };

  void set_not_use_spline(bool flag = true) { not_use_spline = !flag;
  }

  void init_bspline() {
    std::vector<double> tmp(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
      tmp[i] = data[i];
    }

    spline_interp.emplace(tmp, numberLongitudePixels);
    // spline_interp = BSpline2D(tmp, numberLongitudePixels);
  }

  double at_bspline(double lat, double lon) {
    if (not_use_spline) {
      return at(lat, lon);
    }

    if (!spline_interp.has_value()) {
      init_bspline();
    };

    lat *= 3600, lon *= 3600;
    double line = lat - latitudeStartArc, pixel = lon - longitudeStartArc;
    return spline_interp->at(line, pixel);
  };

  void writeAsTiff(const boost::filesystem::path&);

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int) {
    ar& longitudeStartArc& longitudeEndArc& latitudeStartArc& latitudeEndArc&
        numberLongitudePixels& numberLatitudePixels& data;
  };
};

using TSX_indexing = std::vector<std::vector<std::string>>;

// Get TSX tile containing a specific point.
std::string getTileFromTSX(const TSX_indexing& TSX, double lon, double lat);
}  // namespace WHU