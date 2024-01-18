#pragma once

#include <vector>
#include <functional>
#include <array>
#include <fstream>
#include <tuple>

#define _USE_MATH_DEFINES
#include <math.h>

#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/array.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "th_vec.h"
#include "th_serialize.h"
#include "th_state.h"
#include "th_statepolynomial.h"

// Gotcha: Carefully design the header file depencency relationship, and one header per class.
namespace WHU {
    struct Datatake : public ISerializable<Datatake> {
        std::vector<State> states;

        double lambda;

        Timestamp azimuthStartTime;
        Timestamp azimuthEndTime;
        size_t numberOfAzimuthLines;
        double rangeStartDistance;
        double rangeEndDistance;
        size_t numberOfRangePixels;

        int look_direction; // 1 for left_looking, -1 for right_looking

        bool native_doppler;
        std::array<std::array<double, 4>, 3> doppler_coefficients;

        struct NaiveDopplerGetter {
            bool with_value;

            double r1;
            double fdp0, fdp1, fdp2, fdp3;

            double doppler(double range){
                if (!with_value) {
                    return 0;
                };

                double fd;
                double dr, dr2, dr3;

                dr = range - r1;
                dr2 = dr * dr;
                dr3 = dr2 * dr;

                fd = fdp0 + fdp1 * dr + fdp2 * dr2 + fdp3 * dr3;
                return fd;
            }
        };

        NaiveDopplerGetter getNaiveDopplerGetter() const {
            NaiveDopplerGetter tmp;
            if (!native_doppler) {
                tmp.with_value = false;
                return tmp;
            }
            tmp.with_value = true;

            bool with_non_zero_value = false;
            for (int i = 1; i < 3; ++i) {
                for (int j = 0; j < 4; ++j) {
                    with_non_zero_value |= (doppler_coefficients[i][j] != 0);
                }
            }
            assert(with_non_zero_value == false);

            tmp.r1 = (rangeStartDistance + rangeEndDistance) / 2;
            tmp.fdp0 = doppler_coefficients[0][0];
            tmp.fdp1 = doppler_coefficients[0][1];
            tmp.fdp2 = doppler_coefficients[0][2];
            tmp.fdp3 = doppler_coefficients[0][3];
            return tmp;
        };

        friend class boost::serialization::access;

        // TODO: Minimum contructors and no assignment operators are provided. May cause performance issue.
        Datatake(): native_doppler(false) {};
        
        // TODO: too tedious, find a common pattern to reduce codes.
        Datatake(const Datatake& source) {
            states = source.states;
            lambda = source.lambda;
            azimuthStartTime = source.azimuthStartTime;
            azimuthEndTime = source.azimuthEndTime;
            numberOfAzimuthLines = source.numberOfAzimuthLines;
            rangeStartDistance = source.rangeStartDistance;
            rangeEndDistance = source.rangeEndDistance;
            numberOfRangePixels = source.numberOfRangePixels;
            look_direction = source.look_direction;
            for (auto& state : states) {
                state.datatake = this;
            }
            native_doppler = source.native_doppler;
            doppler_coefficients = source.doppler_coefficients;
        };

        double rangePixelSpacing() const {
            return (rangeEndDistance - rangeStartDistance) / (numberOfRangePixels - 1);
        };

        double azimuthLineTime() const {
            return (azimuthEndTime - azimuthStartTime) / (numberOfAzimuthLines - 1);
        };

        BOOST_SERIALIZATION_SPLIT_MEMBER();

        template<class Archive>
        void save(Archive& ar, const unsigned int) const& {
            ar& states& lambda& azimuthStartTime& azimuthEndTime
                & numberOfAzimuthLines& rangeStartDistance& rangeEndDistance& numberOfRangePixels
                & look_direction
                & native_doppler& doppler_coefficients;
        };

        template<class Archive>
        void load(Archive& ar, const unsigned int)& {
            ar& states& lambda& azimuthStartTime& azimuthEndTime
                & numberOfAzimuthLines& rangeStartDistance& rangeEndDistance& numberOfRangePixels
                & look_direction
                & native_doppler& doppler_coefficients;
            for (auto& state : states) {
                state.datatake = this;
            }
        };

        std::tuple<double, double, double, double> getExtent() const;

        WHU::StatePolynomial getSinglePolynomial() const;
    };
};