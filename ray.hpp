#ifndef RAY_HPP
#define RAY_HPP
#include <iostream>
#include "vec3.hpp"


class Ray {
    public:
        Vec3 origin;
        Vec3 direction;
        const double tmin;
        const double tmax;
        Ray() : tmin(0.001), tmax(10000) {};
        Ray(const Vec3& _o, const Vec3& _d) : tmin(0.001), tmax(10000) {
            origin = _o;
            direction = _d;
        };
        ~Ray() {};
};
#endif
