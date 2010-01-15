/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/**
  * This set of classes reads different PNM-files (PGM, PPM, etc)
  *
  * @author André R. Brodtkorb
  */

#ifndef HONEI_UTIL_GUARD_PNM_IMAGE_HH
#define HONEI_UTIL_GUARD_PNM_IMAGE_HH 1

#include <iostream>
#include <exception>
#include <string>

#include <tr1/memory>

namespace honei
{
    namespace util 
    {

        class PNMImageExcpetion : public std::exception {
            public:
                PNMImageExcpetion(const std::string& text_) {
                    text = text_;
                }
                ~PNMImageExcpetion() throw() {}

                virtual const char* what() const throw() {
                    return text.c_str();        
                }
            private:
                std::string text;
        };

        /**
         * Reads an image in pgm-format, and returns values in floating point scaled between 0 and 1.
         *
         * From http://netpbm.sourceforge.net/doc/pgm.html:
         * 
         * Each PGM image consists of the following:
         *
         * 1. A "magic number" for identifying the file type. A pgm image's magic 
         *    number is the two characters "P5".
         * 2. Whitespace (blanks, TABs, CRs, LFs).
         * 3. A width, formatted as ASCII characters in decimal.
         * 4. Whitespace.
         * 5. A height, again in ASCII decimal.
         * 6. Whitespace.
         * 7. The maximum gray value (Maxval), again in ASCII decimal. Must be less 
         *    than 65536, and more than zero.
         * 8. A single whitespace character (usually a newline).
         * 9. A raster of Height rows, in order from top to bottom. Each row 
         *    consists of Width gray values, in order from left to right. Each gray 
         *    value is a number from 0 through Maxval, with 0 being black and Maxval
         *    being white. Each gray value is represented in pure binary by either 1
         *    or 2 bytes. If the Maxval is less than 256, it is 1 byte. Otherwise, 
         *    it is 2 bytes. The most significant byte is first.
         *
         * A row of an image is horizontal. A column is vertical. The pixels in the
         * image are square and contiguous.
         *
         * Each gray value is a number proportional to the intensity of the pixel, 
         * [...]
         * Note that a common variation on the PGM format is to have the gray value
         * be "linear," [...]. pnmgamma takes such a PGM variant as input and produces
         * a true PGM as output.
         *
         * Strings starting with "#" may be comments, the same as with PBM. 
         */
        class PGMImage {
            public:
                PGMImage(size_t width, size_t height);
                ~PGMImage();

                float* getGrayData();
                size_t getWidth();
                size_t getHeight();

                static std::tr1::shared_ptr<PGMImage> read(const char* filename) throw(PNMImageExcpetion);
                static void write(std::tr1::shared_ptr<PGMImage>, const char* filename);

            private:
                PGMImage() {};

            private:
                float* data;
                size_t width;
                size_t height;
        };

        /**
         * Essentially the same as PGMImage, but three channels interleaved. Thus as above except
         * 1. A "magic number" for identifying the file type. A ppm image's magic number is the two
         *    characters "P6". 
         * 9. A raster of Height rows, in order from top to bottom. Each row consists of Width pixels, 
         *    in order from left to right. Each pixel is a triplet of red, green, and blue samples, in 
         *    that order. Each sample is represented in pure binary by either 1 or 2 bytes. If the Maxval
         *    is less than 256, it is 1 byte. Otherwise, it is 2 bytes. The most significant byte is first. 
         */
        class PPMImage {
            public:
                PPMImage(size_t width, size_t height);
                ~PPMImage();

                float* getRedData();
                float* getGreenData();
                float* getBlueData();
                size_t getWidth();
                size_t getHeight();

                static std::tr1::shared_ptr<PPMImage> read(const char* filename) throw(PNMImageExcpetion);
                static void write(std::tr1::shared_ptr<PPMImage>& image, const char* filename) throw(PNMImageExcpetion);

            private:
                PPMImage() {};

            private:
                float* red, * green, * blue;
                size_t width;
                size_t height;
        };






        /**
         * Inline functions
         */
        inline float* PGMImage::getGrayData() { return data; }
        inline size_t PGMImage::getWidth() { return width; }
        inline size_t PGMImage::getHeight() { return height; }

        inline float* PPMImage::getRedData() { return red; }
        inline float* PPMImage::getGreenData() { return green; }
        inline float* PPMImage::getBlueData() { return blue; }
        inline size_t PPMImage::getWidth() { return width; }
        inline size_t PPMImage::getHeight() { return height; }

    }
}
#endif /*PNMIMAGE_H_*/
