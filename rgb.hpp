#ifndef RGB_HPP
#define RGB_HPP
#include <iostream>
#include "vec3.hpp"
inline double clamp(double x, double xmin, double xmax) {
    if(x < xmin) {
        return xmin;
    }
    else if(x > xmax) {
        return xmax;
    }
    else {
        return x;
    }
}
class RGB {
    public:
        double r;
        double g;
        double b;
        RGB() { r = 0; g = 0; b = 0; };
        RGB(double _r) { r = _r; g = _r; b = _r; };
        RGB(double _r, double _g, double _b) { r = _r; g = _g; b = _b; };

        RGB operator+(const RGB& r2) const {
            return RGB(r + r2.r, g + r2.g, b + r2.b);
        };
        RGB operator+(double k) const {
            return RGB(r + k, g + k, b + k);
        };
        RGB operator-(const RGB& r2) const {
            return RGB(r - r2.r, g - r2.g, b - r2.b);
        };
        RGB operator-(double k) const {
            return RGB(r - k, g - k, b - k);
        };
        RGB operator*(const RGB& r2) const {
            return RGB(r * r2.r, g * r2.g, b * r2.b);
        };
        RGB operator*(double k) const {
            return RGB(k * r, k * g, k * b);
        };
        RGB operator/(const RGB& r2) const {
            return RGB(r / r2.r, g / r2.g, b / r2.b);
        };
        RGB operator/(double k) const {
            return RGB(r / k, g / k, b / k);
        };
        RGB operator+=(const RGB& r2) {
            r += r2.r;
            g += r2.g;
            b += r2.b;
        };

        RGB rgb_clamp() const {
            return RGB(clamp(r, 0.0, 1.0), clamp(g, 0.0, 1.0), clamp(b, 0.0, 1.0));
        };

        void print() const {
            std::cout << "(" << r << ", " << g << ", " << b << ")" << std::endl;
        };
};
RGB operator+(double k, const RGB& r2) {
    return RGB(k + r2.r, k + r2.g, k + r2.b);
}
RGB operator-(double k, const RGB& r2) {
    return RGB(k - r2.r, k - r2.g, k - r2.b);
}
RGB operator*(double k, const RGB& r2) {
    return RGB(k * r2.r, k * r2.g, k * r2.b);
}
RGB operator/(double k, const RGB& r2) {
    return RGB(k / r2.r, k / r2.g, k / r2.b);
}
inline RGB pow(const RGB& r, double x) {
    return RGB(std::pow(r.r, x), std::pow(r.g, x), std::pow(r.b, x));
}
#endif
