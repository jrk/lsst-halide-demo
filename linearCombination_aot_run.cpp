//To compile and run:
//Look in linearCombination_aot_compile.cpp

//#define STANDALONE
#ifndef STANDALONE
#include "lsst/afw/image.h"
namespace afwImage = lsst::afw::image;
namespace afwMath  = lsst::afw::math;
#endif

#include <stdio.h>
#include <bitset>
#include "clock.h"
#include <iostream>
using namespace std;

#include "lincombo_aot.h"



void zero_buf_t(buffer_t &buf){
    buf.dev = 0;
    buf.host = 0;
    buf.extent[0] = 0;
    buf.extent[1] = 0;
    buf.extent[2] = 0;
    buf.extent[3] = 0;
    buf.stride[0] = 0;
    buf.stride[1] = 0;
    buf.stride[2] = 0;
    buf.stride[3] = 0;
    buf.min[0] = 0;
    buf.min[1] = 0;
    buf.min[2] = 0;
    buf.min[3] = 0;
    buf.elem_size = 0;
    buf.host_dirty = 0;
    buf.dev_dirty = 0;

}


int main(int argc, char *argv[]) {

    cout << "hello!" << endl;

#ifndef STANDALONE
    auto im = afwImage::MaskedImage<float>("./images/calexp-004207-g3-0123.fits");
    int width = im.getWidth(), height = im.getHeight();
#else
    int width = 2048, height = 1489;
    printf("[no load]");
#endif
    printf("Loaded: %d x %d\n", width, height);

    uint8_t *image = new uint8_t[width*height*4];
    uint8_t *variance = new uint8_t[width*height*4];
    uint8_t *mask = new uint8_t[width*height*2];


#ifndef STANDALONE
    //Read image, converting all three planes to uint8_t arrays
    //for passing to the aot compiled Halide
    float curImage;
    float curVariance;
    uint16_t curMask;

    uint8_t *curImageUInt8Array;
    uint8_t *curVarianceUInt8Array;
    uint8_t *curMaskUInt8Array;

    for (int y = 0; y < height; y++) {
        afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = im.x_at(0, y);
        for (int x = 0; x < width; x++){
            curImage = (*inPtr).image();
            curVariance = (*inPtr).variance();
            curMask = (*inPtr).mask();
            inPtr++;

            curImageUInt8Array = reinterpret_cast<uint8_t*>(&curImage);
            curVarianceUInt8Array = reinterpret_cast<uint8_t*>(&curVariance);
            curMaskUInt8Array = reinterpret_cast<uint8_t*>(&curMask);

            for(int i = 0; i < 4; i++){
                image[(y*width + x)*4 + i] = curImageUInt8Array[i];
                variance[(y*width + x)*4 + i] = curVarianceUInt8Array[i];
            }

            for(int i = 0; i < 2; i++){
                mask[(y*width + x)*2 + i] = curMaskUInt8Array[i];
            }
        }
    }
#endif


    // Have a look in the header file above (it won't exist until you've run
    // lesson_10_generate).

    // It starts with a definition of a buffer_t:
    //
    // typedef struct buffer_t {
    //     uint64_t dev;
    //     uint8_t* host;
    //     int32_t extent[4];
    //     int32_t stride[4];
    //     int32_t min[4];
    //     int32_t elem_size;
    //     bool host_dirty;
    //     bool dev_dirty;
    // } buffer_t;
    //
    // This is how Halide represents input and output images in
    // pre-compiled pipelines. There's a 'host' pointer that points to the
    // start of the image data, some fields that describe how to access
    // pixels, and some fields related to using the GPU that we'll ignore
    // for now (dev, host_dirty, dev_dirty).

    // Let's make some input data to test with:
//    uint8_t input[640 * 480];
//    for (int y = 0; y < 480; y++) {
//        for (int x = 0; x < 640; x++) {
//            input[y * 640 + x] = x ^ (y + 1);
//        }
//    }

    // And the memory where we want to write our output:
    //for every pixel we need 4 image bytes, 4 variance bytes, and 2 mask bytes
    uint8_t image_output[width*height*4];
    uint8_t variance_output[width*height*4];
    uint8_t mask_output[width*height*2];


    // In AOT-compiled mode, Halide doesn't manage this memory for
    // you. You should use whatever image data type makes sense for
    // your application. Halide just needs pointers to it.

    // Now we make a buffer_t to represent our input and output. It's
    // important to zero-initialize them so you don't end up with
    // garbage fields that confuse Halide.
    buffer_t image_buf = {0};
    buffer_t variance_buf = {0};
    buffer_t mask_buf = {0};
    buffer_t image_output_buf = {0};
    buffer_t variance_output_buf = {0};
    buffer_t mask_output_buf = {0};

    // The host pointers point to the start of the image data:
    image_buf.host  = &image[0];
    variance_buf.host  = &variance[0];
    mask_buf.host  = &mask[0];
    image_output_buf.host = &image_output[0];
    variance_output_buf.host = &variance_output[0];
    mask_output_buf.host = &mask_output[0];

    // To access pixel (x, y) in a two-dimensional buffer_t, Halide
    // looks at memory address:

    // host + elem_size * ((x - min[0])*stride[0] + (y - min[1])*stride[1])

    // The stride in a dimension represents the number of elements in
    // memory between adjacent entries in that dimension. We have a
    // grayscale image stored in scanline order, so stride[0] is 1,
    // because pixels that are adjacent in x are next to each other in
    // memory.
    image_buf.stride[0] = variance_buf.stride[0] = mask_buf.stride[0] = 1;
    image_output_buf.stride[0] = variance_output_buf.stride[0] = 1;
    mask_output_buf.stride[0] = 1;

    // stride[1] is the width of the image, because pixels that are
    // adjacent in y are separated by a scanline's worth of pixels in
    // memory.
    image_buf.stride[1] = variance_buf.stride[1] = mask_buf.stride[1] =width;
    image_output_buf.stride[1] = variance_output_buf.stride[1] =width;
    mask_output_buf.stride[1] =width;

    // The extent tells us how large the image is in each dimension.
    image_buf.extent[0] = variance_buf.extent[0] = mask_buf.extent[0] = width;
    image_output_buf.extent[0] = variance_output_buf.extent[0] = width;
    mask_output_buf.extent[0] = width;


    image_buf.extent[1] = variance_buf.extent[1] = mask_buf.extent[1] = height;
    image_output_buf.extent[1] = variance_output_buf.extent[1] = height;
    mask_output_buf.extent[1] = height;


    // We'll leave the mins as zero. This is what they typically
    // are. The host pointer points to the memory location of the min
    // coordinate (not the origin!).  See lesson 6 for more detail
    // about the mins.

    // The elem_size field tells us how many bytes each element
    // uses. For the 8-bit image we use in this test it's one.
    image_buf.elem_size = variance_buf.elem_size = 4;
    mask_buf.elem_size = 2;

    image_output_buf.elem_size = variance_output_buf.elem_size = 4;
    mask_output_buf.elem_size = 2;

    // To avoid repeating all the boilerplate above, We recommend you
    // make a helper function that populates a buffer_t given whatever
    // image type you're using.

    // Now that we've setup our input and output buffers, we can call
    // our function. Looking in the header file, it's signature is:

    // int test_aot(buffer_t *_input, const int32_t _offset, buffer_t *_brighter);

    // The return value is an error code. It's zero on success.

    int error = lincombo_aot(&image_buf, &variance_buf, &mask_buf,
                    &image_output_buf, &variance_output_buf, &mask_output_buf);

    if (error) {
        printf("Halide returned an error: %d\n", error);
        return -1;
    }


    // Benchmark the pipeline.
    double mean = 0;
    double min;
    double max;
    int numberOfRuns = 5;
    for (int i = 0; i < numberOfRuns; i++) {
        double t1 = current_time();
        error = lincombo_aot(&image_buf, &variance_buf, &mask_buf,
                    &image_output_buf, &variance_output_buf, &mask_output_buf);
        double t2 = current_time();
        double curTime = (t2-t1);
        mean += curTime;
        if(i == 0){
            min = curTime;
            max = curTime;
        }
        else{
            if(curTime < min)
                min = curTime;
            if(curTime > max)
                max = curTime;
        }
    }
    mean = mean/numberOfRuns;

    std::cout << "Mean Time: " << mean << ", Min = " <<
    min << ", Max = " << max << ", with " << numberOfRuns <<
    " runs" << '\n';

#ifndef STANDALONE
    bool writePlanesSeparately = true;
    if(!writePlanesSeparately){
        //write image out
        auto imOut = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel>(im.getWidth(), im.getHeight());
        for (int y = 0; y < imOut.getHeight(); y++) {
        	afwImage::MaskedImage<float, lsst::afw::image::MaskPixel,
                lsst::afw::image::VariancePixel>::x_iterator inPtr = imOut.x_at(0, y);
            for (int x = 0; x < imOut.getWidth(); x++){
                curImageUInt8Array = image_output + 4*(y*width + x);
                curVarianceUInt8Array = variance_output + 4*(y*width + x);
                curMaskUInt8Array = mask_output + 2*(y*width + x);

                curImage = *(reinterpret_cast<float*>(curImageUInt8Array));
                curVariance = *(reinterpret_cast<float*>(curVarianceUInt8Array));
                curMask = *(reinterpret_cast<uint16_t*>(curMaskUInt8Array));

                afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel,
                    lsst::afw::image::VariancePixel> curPixel(curImage, curMask, curVariance);

            	(*inPtr) = curPixel;
            	inPtr++;

            }
        }

    	imOut.writeFits("./halideCleanLinearCombination5x5.fits");
    }
    else{
        //write three planes separately
        auto imOut = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(im.getWidth(), im.getHeight());
        for (int y = 0; y < imOut.getHeight(); y++) {
            afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = imOut.x_at(0, y);

            for (int x = 0; x < imOut.getWidth(); x++){
                curImageUInt8Array = image_output + 4*(y*width + x);
                curImage = *(reinterpret_cast<float*>(curImageUInt8Array));
                afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel,
                lsst::afw::image::VariancePixel> curPixel(curImage, 0, 0);
                (*inPtr) = curPixel;
                inPtr++;

            }
        }

        auto varOut = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(im.getWidth(), im.getHeight());
        for (int y = 0; y < imOut.getHeight(); y++) {
            afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = varOut.x_at(0, y);

            for (int x = 0; x < imOut.getWidth(); x++){
                curVarianceUInt8Array = variance_output + 4*(y*width + x);
                curVariance = *(reinterpret_cast<float*>(curVarianceUInt8Array));

                afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel,
                lsst::afw::image::VariancePixel> curPixel(curVariance, 0, 0);
                (*inPtr) = curPixel;
                inPtr++;

            }
        }

        auto maskOutPlane = afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(im.getWidth(), im.getHeight());
        for (int y = 0; y < imOut.getHeight(); y++) {
            afwImage::MaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::x_iterator inPtr = maskOutPlane.x_at(0, y);

            for (int x = 0; x < imOut.getWidth(); x++){
                curMaskUInt8Array = mask_output + 2*(y*width + x);
                curMask = *(reinterpret_cast<uint16_t*>(curMaskUInt8Array));

                afwImage::pixel::SinglePixel<float, lsst::afw::image::MaskPixel,
                lsst::afw::image::VariancePixel> curPixel(curMask, 0, 0);
                (*inPtr) = curPixel;
                inPtr++;

            }
        }

        imOut.writeFits("./halideLinComboImage5x5.fits");
        varOut.writeFits("./halideLinComboVar5x5.fits");
        maskOutPlane.writeFits("./halideLinComboMask5x5.fits");
    }
#endif
}


