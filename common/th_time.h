#pragma once

#include <string>
#include <cassert>

#include <fmt/format.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace WHU {
    struct Timestamp {
        int year;
        int month;
        int day;
        double seconds;

        static Timestamp from_TSX_timestamp_string(const std::string& timestamp) {
            Timestamp tmp = {
                std::stoi(timestamp.substr(0, 4)),
                std::stoi(timestamp.substr(5, 2)),
                std::stoi(timestamp.substr(8, 2)),
                std::stoi(timestamp.substr(11, 2)) * 3600 +
                std::stoi(timestamp.substr(14, 2)) * 60 +
                std::stod(timestamp.substr(17, 9))
            };

            return tmp;
        };

        void assert_same_day(const Timestamp& other) const {
            assert(year == other.year && month == other.month && day == other.day);
        }

        double operator-(const Timestamp& other) const {
            assert_same_day(other);
            return seconds - other.seconds;
        };

        Timestamp operator+(double _seconds) const {
            assert(seconds + _seconds< 86400);
            return Timestamp(year, month, day, seconds + _seconds);
        }

        Timestamp& operator+=(double _seconds) {
            assert(seconds + _seconds < 86400);
            seconds += _seconds;
            return *this;
        }

        operator double() const {
            return seconds;
        } 

        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive& ar, const unsigned int) {
            ar& year& month& day& seconds;
        };
    };
};