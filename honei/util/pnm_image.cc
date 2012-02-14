/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <honei/util/pnm_image.hh>
#include <honei/util/stringify.hh>

#include <iostream>
#include <fstream>
#include <inttypes.h>
//#include <boost/cstdint.hpp>

using namespace honei;
using namespace util;
using namespace std;
//using namespace boost;

inline void stripComment(ifstream& file) {
    char c;
    file >> c;
    if (c == '#')
        while (file.get() != '\n');
    //file.ignore(65336, '\n');
    else
        file.putback(c);
}

/**
 * Reads the header of a PNM-image
 * @templateparam magic Expected magic number.
 * @param file The (binary) ifstream to read from.
 * @param bpc number of *byte* per pixel channel
 * @param max_value maximum value in image
 * @param width Image width in pixels
 * @param height Image height in pixels
 */
    template <char magic>
inline void readPNMHeader(ifstream& file,
        unsigned int& bpc, unsigned int& max_value,
        unsigned int& width, unsigned int& height)
throw(PNMImageExcpetion) {
    char m[2];

    //1. "magic number"
    stripComment(file);
    file >> m[0];//P
    file >> m[1];//1..6
    if (m[0]!= 'P' || m[1] != magic)
        throw PNMImageExcpetion("Wrong magic number");

    //3-5. Width and height
    stripComment(file);
    file >> width;
    stripComment(file);
    file >> height;

    //7 maximum gray value
    stripComment(file);
    file >> max_value;
    if (max_value > 65535)
        throw PNMImageExcpetion("Wrong max value");
    bpc = (max_value < 256) ? 1 : 2;

    //8 skip whitespace
    stripComment(file);
}

/**
 * Reads the header of a PNM-image
 * @templateparam magic Expected magic number.
 * @param file The (binary) ifstream to read from.
 * @param bpc number of *byte* per pixel channel
 * @param max_value maximum value in image
 * @param width Image width in pixels
 * @param height Image height in pixels
 */
    template <char magic>
inline void writePNMHeader(ofstream& file,
        unsigned int& bpc, unsigned int& max_value,
        unsigned int& width, unsigned int& height)
throw(PNMImageExcpetion) {
    max_value = 250;//65535;
    bpc = 1;

    //1. "magic number"
    file << 'P' << magic << endl;

    //comment
    file << "#Written by babrodtks really cool PNMwriter" << endl;

    //3-5. Width and height
    file << width << " " << height << endl;

    //7 maximum gray value
    file << max_value << endl;
}


inline uint16_t read16BitValue(uint8_t* data) {
    uint16_t pixel;
    pixel  = data[1] << 8;
    pixel |= data[0];
    return pixel;
}

inline void fillBuffer(ifstream& file, unsigned int size, uint8_t*& buffer) {
    unsigned int read_bytes = 0;
    unsigned int i = 0;
    while (read_bytes != size && ++i<100) {
        file.read((char*) buffer, size);
        read_bytes += file.gcount();
    }

    if (read_bytes != size)
        throw PNMImageExcpetion("Unable to read pixel data properly");
}

inline void writeBuffer(ofstream& file, unsigned int size, uint8_t*& buffer) {
    file.write((char*) buffer, size);

    if (!file.good())
        throw PNMImageExcpetion("Unable to write pixel data properly");
}

inline void openIfstream(const char*& filename, ifstream& file) {
    file.open(filename, ios::in | ios::out | ios::binary);
    if (!file)
        throw PNMImageExcpetion("Unable to open file " + stringify(filename));
}

inline void openOfstream(const char*& filename, ofstream& file) {
    file.open(filename, ios::out | ios::binary);
    if (!file)
        throw PNMImageExcpetion("Unable to open file " + stringify(filename));
}

PGMImage::PGMImage(size_t width, size_t height) {
    data = new float[width*height];
    this->width = width;
    this->height = height;
}

PGMImage::~PGMImage() {
    delete [] data;
}

std::tr1::shared_ptr<PGMImage> PGMImage::read(const char* filename) throw(PNMImageExcpetion) {
    std::tr1::shared_ptr<PGMImage> img;

    ifstream file;
    float* data;
    uint8_t* buffer;

    unsigned int width;
    unsigned int height;
    unsigned int max_value;
    unsigned int bpc;

    //Open file and read header
    openIfstream(filename, file);
    readPNMHeader<'5'>(file, bpc, max_value, width, height);

    //Allocate data and read into buffer
    img.reset(new PGMImage(width, height));
    buffer = new uint8_t[width*height*bpc];
    fillBuffer(file, width*height*bpc, buffer);

    //Convert to float representation between 0 and 1
    data = img->getGrayData();
    if (bpc == 2) {
        for(unsigned long i(0) ; i < width * height ; ++i) {
            data[i] = read16BitValue(&buffer[2*i]) / (float) max_value;
        }
    }
    else {
        for(unsigned long i(0) ; i < width * height ; ++i) {
            data[i] = buffer[i] / (float) max_value;
        }
    }

    delete [] buffer;
    file.close();

    return img;
}


PPMImage::PPMImage(size_t width, size_t height) {
    red = new float[width*height];
    green = new float[width*height];
    blue = new float[width*height];
    this->width = width;
    this->height = height;
}

PPMImage::~PPMImage() {
    delete [] red;
    delete [] green;
    delete [] blue;
}

std::tr1::shared_ptr<PPMImage> PPMImage::read(const char* filename) throw(PNMImageExcpetion) {
    std::tr1::shared_ptr<PPMImage> img;

    ifstream file;
    float* data[3];
    uint8_t* buffer;

    unsigned int width;
    unsigned int height;
    unsigned int max_value;
    unsigned int bpc;

    //Open file and read header
    openIfstream(filename, file);
    readPNMHeader<'6'>(file, bpc, max_value, width, height);

    //Allocate data and read into buffer
    img.reset(new PPMImage(width, height));
    buffer = new uint8_t[3*width*height*bpc];
    fillBuffer(file, 3*width*height*bpc, buffer);

    //Convert to float representation between 0 and 1
    data[0] = img->getRedData();
    data[1] = img->getGreenData();
    data[2] = img->getBlueData();
    if (bpc == 2) {
        for(unsigned long i(0) ; i < width * height ; ++i) {
            for (int j=0; j<3; ++j) {
                data[j][i] = read16BitValue(&buffer[6*i+2*j]) / (float) max_value;
            }
        }
    }
    else {
        for(unsigned long i(0) ; i < width * height ; ++i) {
            for (int j=0; j<3; ++j) {
                data[j][i] = buffer[3*i+j] / (float) max_value;
            }
        }
    }

    delete [] buffer;
    file.close();

    return img;
}


void PPMImage::write(std::tr1::shared_ptr<PPMImage>& img, const char* filename) throw(PNMImageExcpetion) {
    ofstream file;
    float* data[3];
    uint8_t* buffer;

    unsigned int width = img->getWidth();
    unsigned int height = img->getHeight();
    unsigned int max_value;
    unsigned int bpc;

    //Open file and read header
    openOfstream(filename, file);
    writePNMHeader<'6'>(file, bpc, max_value, width, height);

    //Allocate data and read into buffer
    buffer = new uint8_t[3*width*height*bpc];

    //Convert to float representation between 0 and 1
    data[0] = img->getRedData();
    data[1] = img->getGreenData();
    data[2] = img->getBlueData();
    if (bpc == 2) {
        for(unsigned long i(0); i < width * height ; ++i) {
            for (int j=0; j<3; ++j) {
                buffer[6*i  ] = data[0][i]*max_value;
                buffer[6*i+1] = data[0][i+1]*max_value;
                buffer[6*i+2] = data[1][i]*max_value;
                buffer[6*i+3] = data[1][i+1]*max_value;
                buffer[6*i+4] = data[2][i]*max_value;
                buffer[6*i+5] = data[2][i+1]*max_value;
            }
        }
    }
    else {
        for(unsigned long i(0) ; i < width * height ; ++i) {
            for (int j=0; j<3; ++j) {
                buffer[3*i  ] = data[0][i]*max_value;
                buffer[3*i+1] = data[1][i]*max_value;
                buffer[3*i+2] = data[2][i]*max_value;
            }
        }
    }

    writeBuffer(file, 3*width*height*bpc, buffer);

    delete [] buffer;

    file.close();
}
