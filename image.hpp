#ifndef IMAGE_HPP
#define IMAGE_HPP
#include <iostream>
#include <fstream>
#include <string>
#include "rgb.hpp"
class Image {
    public:
        RGB* color;
        int width;
        int height;
        Image(int _width, int _height) {
            width = _width;
            height = _height;
            color = new RGB[width*height];
        };
        ~Image() {
            delete[] color;
        };
        void setPixel(int i, int j, RGB col) { color[i + j*width] = col; };
        RGB getPixel(int i, int j) { return color[i + j*width]; };
        void ppm_output(std::string filename) {
            std::ofstream file(filename);

            file << "P3\n";
            file << width << " " << height << "\n";
            file << 255 << "\n";

            for(int j = 0; j < height; j++) {
                for(int i = 0; i < width; i++) {
                    RGB col = color[i + j*width];
                    col.rgb_clamp();
                    int r = clamp((int)(255*col.r), 0, 255);
                    int g = clamp((int)(255*col.g), 0, 255);
                    int b = clamp((int)(255*col.b), 0, 255);
                    file << r << " " << g << " " << b << "\n";
                }
            }
            file.close();
        };
        void divide(double k) {
            for(int i = 0; i < width; i++)
                for(int j = 0; j < height; j++)
                    color[i + j*width] = color[i + j*width]/k;
        };
        void gamma_correlation() {
            for(int i = 0; i < width; i++)
                for(int j = 0; j < height; j++)
                    color[i + j*width] = pow(color[i + j*width], 1.0/2.2);
        };
};
#endif
