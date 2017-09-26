#include <iostream>
#include <string>
#include <cmath>
#include <omp.h>

#include "vec3.hpp"
#include "ray.hpp"
#include "camera.hpp"
#include "rgb.hpp"
#include "image.hpp"


const int W = 800;
const int H = 800;
const Vec3 sunDir = Vec3(0, 1, -1).normalize();


inline double rayleigh_phase(const Vec3& rayDir, const Vec3& lightDir) {
    double mu = rayDir.dot(lightDir);
    return 3.0/(16.0*M_PI) * (1.0 + mu*mu);
}
inline double mie_phase(const Vec3& rayDir, const Vec3& lightDir, double g) {
    double mu = rayDir.dot(lightDir);
    return 3.0/(8.0*M_PI) * ((1.0 - g*g)*(1.0 + mu*mu))/((2.0 + g*g)*std::pow((1.0 + g*g -2.0*g*mu), 1.5));
}


struct Hit {
    double t;
    Vec3 hitPos;
    Vec3 hitNormal;
    Hit() {};
};
bool raySphereIntersection(const Ray& ray, const Vec3& center, double radius, Hit& res) {
    double b = ray.direction.dot(ray.origin - center);
    double c = (ray.origin - center).length2() - radius*radius;
    double D = b*b - c;
    if(D < 0) {
        return false;
    }
    double t;
    double t1 = -b + std::sqrt(D);
    double t2 = -b - std::sqrt(D);
    if(t1 > 0 && t2 > 0) {
        t = std::min(t1, t2);
    }
    else if(t1 > 0 && t2 < 0) {
        t = t1;
    }
    else if(t1 < 0 && t2 > 0) {
        t = t2;
    }
    else {
        return false;
    }
    res.t = t;
    res.hitPos = ray.origin + t*ray.direction;
    res.hitNormal = (res.hitPos - center).normalize();
    return true;
}


const Vec3 earth_origin = Vec3(0);
const double atmos_radius = 6420 * 1000;
const double earth_radius = 6360 * 1000;
const double rayleigh_scaleheight = 7400;
const double mie_scaleheight = 1200;
const double sun_intensity = 20;
const RGB beta_rayleigh = RGB(3.8e-6, 13.5e-6, 33.1e-6);
const RGB beta_mie = RGB(21e-6);
const int view_samples = 20;
const int light_samples = 10;
RGB sky(const Ray& ray) {
    Hit res;
    if(!raySphereIntersection(ray, earth_origin, atmos_radius, res)) {
        std::cout << "View ray didn't hit the outer bound of atmosphere" << std::endl;
        return RGB(0);
    }
    double atmos_dist = res.t;
    
    double ds = atmos_dist/(double)view_samples;
    double phase_rayleigh = rayleigh_phase(ray.direction, sunDir);
    double phase_mie = mie_phase(ray.direction, sunDir, 0.74);
    double optical_depth_rayleigh = 0.0;
    double optical_depth_mie = 0.0;
    RGB color_rayleigh = RGB(0);
    RGB color_mie = RGB(0);
    for(int i = 0; i < view_samples; i++) {
        double t = ds * (i + 0.5);
        Vec3 samplePos = ray.origin + t*ray.direction;
        double sea_level = (samplePos - earth_origin).length() - earth_radius;
        double h_rayleigh = std::exp(-sea_level/rayleigh_scaleheight) * ds;
        double h_mie = std::exp(-sea_level/mie_scaleheight) * ds;
        optical_depth_rayleigh += h_rayleigh;
        optical_depth_mie += h_mie;

        //Direct Illumination Sampling
        Ray lightRay = Ray(samplePos, sunDir);
        Hit res2;
        if(!raySphereIntersection(lightRay, earth_origin, atmos_radius, res2)) {
            std::cout << "Light ray didn't hit the outer bound of atmosphere" << std::endl;
            return RGB(0);
        }
        double atmos_dist_lightray = res2.t;
        double ds_light = atmos_dist_lightray/(double)light_samples;
        double optical_depth_rayleigh_light = 0.0;
        double optical_depth_mie_light = 0.0;
        for(int j = 0; j < light_samples; j++) {
            double t2 = ds_light * (j + 0.5);
            Vec3 samplePos_light = lightRay.origin + t2*lightRay.direction;
            double sea_level_light = (samplePos_light - earth_origin).length() - earth_radius;
            if(sea_level_light < 0) {
                //std::cout << "light ray penetrated the earth" << std::endl;
                break;
            }
            optical_depth_rayleigh_light += std::exp(-sea_level_light/rayleigh_scaleheight) * ds_light;
            optical_depth_mie_light += std::exp(-sea_level_light/mie_scaleheight) * ds_light;

            if(j == light_samples - 1) {
                RGB tau = beta_rayleigh * (optical_depth_rayleigh + optical_depth_rayleigh_light) + beta_mie * 2.1 * (optical_depth_mie + optical_depth_mie_light);
                RGB attenuation(std::exp(-tau.r), std::exp(-tau.g), std::exp(-tau.b));
                color_rayleigh += attenuation * h_rayleigh;
                color_mie += attenuation * h_mie;
            }
        }
    }
    return (color_rayleigh*beta_rayleigh*phase_rayleigh + color_mie*beta_mie*phase_mie) * sun_intensity;
}


std::string percentage(double x, double max) {
    return std::to_string(x/max*100) + "%";
}
std::string progressbar(double x, double max) {
    int max_count = 40;
    int cur_count = (int)(x/max*max_count);
    std::string str;
    str += "[";
    for(int i = 0; i < cur_count; i++) {
        str += "#";
    }
    for(int i = 0; i < (max_count - cur_count - 1); i++) {
        str += " ";
    }
    str += "]";
    return str;
}


int main() {
    Image img(W, H);

    Vec3 camPos = Vec3(0, earth_radius + 1, 0);
    Vec3 camFront = Vec3(0, 0, -1).normalize();
    PinholeCamera cam(camPos, camFront, 90.0);

    #pragma omp parallel for schedule(dynamic, 1)
    for(int i = 0; i < W; i++) {
        for(int j = 0; j < H; j++) {
            double u = (2.0*i - W)/H;
            double v = (2.0*j - H)/H;
            Ray ray = cam.getRay(u, v);
            RGB col = sky(ray);
            img.setPixel(i, j, col);
            if(omp_get_thread_num() == 0) {
                double progress = (double)(j + H*i)/(W*H);
                std::cout << progressbar(progress, 1.0) << percentage(progress, 1.0) << "\r" << std::flush;
            }
        }
    }
    img.gamma_correlation();
    img.ppm_output("output.ppm");
}
