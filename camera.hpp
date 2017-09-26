#ifndef CAMERA_HPP
#define CAMERA_HPP

#include <cmath>
#include "vec3.hpp"
#include "ray.hpp"

class PinholeCamera {
    public:
        Vec3 camPos;
        Vec3 camFront;
        Vec3 camRight;
        Vec3 camUp;
        double d;
        PinholeCamera() {};
        PinholeCamera(const Vec3& _camPos, const Vec3& _camFront, double fov) {
            camPos = _camPos;
            camFront = _camFront;
            camRight = camFront.cross(Vec3(0, 1, 0));
            camUp = -camRight.cross(camFront);
            camFront = camFront.normalize();
            camRight = camRight.normalize();
            camUp = camUp.normalize();
            d = 1.0/tan(fov/2.0);
        }
        Ray getRay(double u, double v) const {
            Vec3 rayDir = (d*camFront + u*camRight + v*camUp).normalize();
            return Ray(camPos, rayDir);
        };
        Ray getRay2(double u, double v) const {
            if(u*u + v*v > 1) {
                return Ray(Vec3(0), Vec3(0));
            };
            double y = std::sqrt(1.0 - u*u - v*v);
            Vec3 rayDir = Vec3(u, y, v).normalize();
            return Ray(camPos, rayDir);
        };
};


#endif
