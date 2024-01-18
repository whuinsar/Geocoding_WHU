#pragma once

#include <tuple>

#include "th_datatake.h"

namespace WHU {
GeodeticCoordinate RDGeolocation(const State& state,
                                 double dopplerFrequency,
                                 double rho,
                                 double height,
                                 double thres = 0.0001);

GeodeticCoordinate cartesian2geodetic(const Position& p);

Position geodetic2cartesian(const GeodeticCoordinate& p);

void set_longitude_from_zero(bool flag = true);
}  // namespace WHU