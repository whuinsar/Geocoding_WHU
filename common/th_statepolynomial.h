#pragma once

#include "th_interp.h"
#include "th_state.h"
#include "th_time.h"

namespace WHU {
    struct Datatake;
    
    struct StatePolynomial : PolynomialNT<6> {
        Datatake* datatake;

        static StatePolynomial fromObjects(const std::vector<State> states, size_t start, size_t count){
            return StatePolynomial(states, start, count);
        };
        
        StatePolynomial(const std::vector<State> states, size_t start, size_t count) :
            PolynomialNT<6>(PolynomialNT<6>::fromObjects<State>(states, start, count)),
            datatake(states[0].datatake) {};

        // Too ugly, not prudent enough when designing the interface.
        State at(const Timestamp& t) const {
            auto tmp = PolynomialNT<6>::at(t.seconds);
            return State{
                {tmp[0], tmp[1], tmp[2]},
                {tmp[3], tmp[4], tmp[5]},
                t,
                datatake
            };
        }
    };

}