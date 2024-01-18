#pragma once

namespace WHU {
    // Always use WGS-84 ellipsoid parameters.
    static inline double wgs_ra = 6378137,
        wgs_e2 = 6.69437999014e-3,
        ra = wgs_ra,
        e2 = wgs_e2;

    static inline double speedOfLight = 299792458;

    static inline double rad2deg = 57.2957795131;
    static inline double deg2rad = 0.0174532925199;

    static inline double deg2meter = 111319.490793273;
}
