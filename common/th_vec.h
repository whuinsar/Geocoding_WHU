#pragma once

#include <cmath>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "th_time.h"

namespace WHU {
    struct Vec3 {
        double x;
        double y;
        double z;

        double dot(const Vec3& other) const {
            return x * other.x + y * other.y + z * other.z;
        }

        double norm() const{
            return std::sqrt(dot(*this));
        }

        Vec3 operator/(double k) const {
            assert(k != 0);
            return Vec3{
                x / k,
                y / k,
                z / k
            };
        }

        Vec3 operator*(double k) const {
            return Vec3{
                x * k,
                y * k,
                z * k
            };
        }

        Vec3 operator-(const Vec3& other) const {
            return Vec3{
                x - other.x,
                y - other.y,
                z - other.z
            };
        }

        double distance2(const Vec3& other) const{
            return ((*this) - other).norm();
        }

        Vec3 cross(const Vec3& other) const {
            return Vec3{
                y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x
            };
        };
        
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive& ar, const unsigned int) {
            ar& x& y& z;
        };
    };

    struct GeodeticCoordinate {
        double lon;
        double lat;
        double height;
    };

    using Position = Vec3;
    using Velocity = Vec3;
};