#pragma once

#include "th_vec.h"

namespace WHU {
    struct Datatake;

    struct State {
        Position pos;
        Velocity vel;
        Timestamp time;
        Datatake* datatake;

        static auto getTimeAccessor() {
            auto tmp = [](const State& state) {
                return state.time;
            };
            return tmp;
        };

        static auto getDataAccessors() {
            auto get_pos_x = [](const State& state) { return state.pos.x; };
            auto get_pos_y = [](const State& state) { return state.pos.y; };
            auto get_pos_z = [](const State& state) { return state.pos.z; };
            auto get_vel_x = [](const State& state) { return state.vel.x; };
            auto get_vel_y = [](const State& state) { return state.vel.y; };
            auto get_vel_z = [](const State& state) { return state.vel.z; };
            std::array<std::function<double(const State&)>, 6> tmp = {
                get_pos_x, get_pos_y, get_pos_z,
                get_vel_x, get_vel_y, get_vel_z
            };
            return tmp;
        };

        double estimateTimeDelta(const Position& p) {
            auto look_vector = p - pos;
            auto lv_at_vel = vel.dot(look_vector) / vel.dot(vel);
            return lv_at_vel;
        };

        double estimateTimeDeltaWithDoppler(const Position& p, double doppler);

        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive& ar, const unsigned int) {
            ar& pos& vel& time;
        };
    };
};