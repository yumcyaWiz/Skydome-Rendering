#ifndef VEC3_HPP
#define VEC3_HPP
#include <iostream>
#include <cmath>


class Vec3 {
    public:
        double x;
        double y;
        double z;

        Vec3() { x = 0; y = 0; z = 0; };
        Vec3(double _x) { x = _x; y = _x; z = _x; };
        Vec3(double _x, double _y, double _z) { x = _x; y = _y; z = _z; };
        ~Vec3() {};

        Vec3 operator+(const Vec3& v2) const {
            return Vec3(x + v2.x, y + v2.y, z + v2.z);
        };
        Vec3 operator+(double k) const {
            return Vec3(x + k, y + k, z + k);
        };
        Vec3 operator-(const Vec3& v2) const {
            return Vec3(x - v2.x, y - v2.y, z - v2.z);
        };
        Vec3 operator-(double k) const {
            return Vec3(x - k, y - k, z - k);
        };
        Vec3 operator-() const {
            return Vec3(-x, -y, -z);
        }
        Vec3 operator*(const Vec3& v2) const {
            return Vec3(x * v2.x, y * v2.y, z * v2.z);
        };
        Vec3 operator*(double k) const {
            return Vec3(k * x, k * y, k * z);
        };
        Vec3 operator/(const Vec3& v2) const {
            return Vec3(x / v2.x, y / v2.y, z / v2.z);
        };
        Vec3 operator/(double k) const {
            return Vec3(x / k, y / k, z / k);
        };

        double length() const {
            return std::sqrt(x*x + y*y + z*z);
        };
        double length2() const {
            return x*x + y*y + z*z;
        };

        double dot(const Vec3& v2) const {
            return x*v2.x + y*v2.y + z*v2.z;
        };
        Vec3 cross(const Vec3& v2) const {
            return Vec3(y*v2.z - z*v2.y, z*v2.x - x*v2.z, x*v2.y - y*v2.x);
        };
        
        Vec3 normalize() const {
            return Vec3(x, y, z)/this->length();
        }

        void print() const {
            std::cout << "(" << x << ", " << y << ", " << z << ")" << std::endl;
        }
};
inline Vec3 operator+(double k, const Vec3& v2) {
    return Vec3(k + v2.x, k + v2.y, k + v2.z);
};
inline Vec3 operator-(double k, const Vec3& v2) {
    return Vec3(k - v2.x, k - v2.y, k - v2.z);
};
inline Vec3 operator*(double k, const Vec3& v2) {
    return Vec3(k * v2.x, k * v2.y, k * v2.z);
};
inline Vec3 operator/(double k, const Vec3& v2) {
    return Vec3(k / v2.x, k / v2.y, k / v2.z);
};
#endif
