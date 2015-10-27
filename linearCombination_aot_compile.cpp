/*
 *
 * This file serves of an example of how to perform "ahead of time"
 * (AOT) compilation of Halide pipelines for later use by other
 * applications.
 *
 * The code defines a Halide pipeline and its accompanying schedule,
 * and then compiles the resulting pipeline to a object file and
 * associated C header file. The object file can then be linked into
 * any application that wishes to use the resulting Halide pipeline.
 *
 * The Halide pipeline uses a spatial varying convolution kernel
 * formed from the weighted combination of 5 basis kernels.
 *
 */

#include <stdio.h>
#include <Halide.h>
#include <bitset>
using namespace std;
using namespace Halide;

int main(int argc, char *argv[]) {

    const int num_kernels = 5;

    ImageParam image(type_of<float>(), 2);
    ImageParam variance(type_of<float>(), 2);
    ImageParam mask(type_of<uint16_t>(), 2);

    ImageParam polynomialCoefficients(type_of<float>(), 2); //array with dimension [10][5]
    //Kernel parameters array with dimensions [3][5]
    //first dimension: sigmaX, sigmaY, theta
    ImageParam kerParams(type_of<float>(), 2);


    /*
     *
     * First define the Halide pipeline.  This is the functional
     * specification of the image processing algorithm.
     *
     */

    //Kernel has dimensions (boundingBox*2 + 1) x (boundingBox*2 + 1)
    int boundingBox = 2;
    Var x, y, i, j, y0, yi;
    float pi = 3.14159265359f;


    Func polynomials[num_kernels];
    for(int k = 0; k < num_kernels; k++){
        polynomials[k](x, y) = polynomialCoefficients(0,k) + 
            polynomialCoefficients(1,k)*x + polynomialCoefficients(2,k)*y +
            polynomialCoefficients(3,k)*x*x + polynomialCoefficients(4,k)*x*y + 
            polynomialCoefficients(5,k)*y*y + polynomialCoefficients(6,k)*x*x*x +
            polynomialCoefficients(7,k)*x*x*y + polynomialCoefficients(8,k)*x*y*y
            + polynomialCoefficients(9,k)*y*y*y;
    }

    Func kernels[num_kernels];
    for(int k = 0; k < num_kernels; k++){
        kernels[k](i, j) = exp(-((i*cos(kerParams(2,k)) + j*sin(kerParams(2,k)))*
                    (i*cos(kerParams(2,k)) + j*sin(kerParams(2,k))))
                    / (2*kerParams(0,k)*kerParams(0,k))
                    -((j*cos(kerParams(2,k)) - i*sin(kerParams(2,k)))*
                    (j*cos(kerParams(2,k)) - i*sin(kerParams(2,k))))
                    / (2*kerParams(1,k)*kerParams(1,k))) /
                    (2.0f*pi*kerParams(0,k)*kerParams(1,k));
    }


    //Compute output image plane
    Func image_bounded ("image_bounded");
    image_bounded = BoundaryConditions::repeat_edge(image);
    Expr blur_image_help = 0.0f;
    Expr norm = 0.0f;

    //Compute output variance plane
    Func variance_bounded ("variance_bounded");
    variance_bounded = BoundaryConditions::repeat_edge(variance);
    Func blurVariance ("blurVariance");
    Expr blur_variance_help = 0.0f;

    //Compute output mask plane
    Func mask_bounded ("mask_bounded");
    mask_bounded = BoundaryConditions::repeat_edge(mask);
    Func maskOut ("maskOut");
    Expr maskOutHelp = cast<uint16_t>(0);

    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
            Expr curKernelVal = 0.0f;
            for(int k = 0; k < num_kernels; k++){
                curKernelVal += polynomials[k](x, y)*kernels[k](i, j);
            }
            blur_image_help += image_bounded(x + i, y + j)*curKernelVal;
            blur_variance_help += variance_bounded(x + i, y + j)*curKernelVal*curKernelVal;
            maskOutHelp = select(curKernelVal == 0.0f, maskOutHelp,
                                maskOutHelp | mask_bounded(x + i, y + j));
            norm += curKernelVal;
        }
    }
    blur_image_help = blur_image_help/norm;
    blur_variance_help = blur_variance_help/(norm*norm);

    //set the image edges
    //image edge should be NAN, but this oddly isn't working 
    Expr setEdge = x < boundingBox || y < boundingBox ||
                   x > (image.width() - 1 - boundingBox) ||
                   y > (image.height() - 1 - boundingBox);
    blur_image_help = select(setEdge, INFINITY, blur_image_help); 
    blur_variance_help = select(setEdge, INFINITY, blur_variance_help);
    maskOutHelp = select(setEdge, 16, maskOutHelp);

    //Evaluate image, mask, and variance planes concurrently using a tuple
    Func combined_output ("combined_output");
    combined_output(x, y) = Tuple(blur_image_help, blur_variance_help, maskOutHelp);

    /*
     *
     * Here is the definition of the Halide schedule.  The schedule
     * describes how to implement the kernel described above.
     *
     */

    // Split the y coordinate of the output into strips of 32 scanlines:
    combined_output.split(y, y0, yi, 32);

    // Compute the strips in parallel using all the machines
    // cores. (You can think of processing a strip as a job that
    // placed task queue, and serviced by a bunch of worker threads).
    combined_output.parallel(y0);

    // Vectorize across x loop by a factor of eight.
    combined_output.vectorize(x, 8);


    /*
     *
     * The statement below generates the .o file (and .h) for the
     * Halide kernel that can be linked against by other applications.
     *
     */


    std::vector<Argument> args = {image, variance, mask, polynomialCoefficients, kerParams};
    combined_output.compile_to_file("lincombo_aot", args);

    printf("Halide pipeline has been compiled.  You can now link against the resulting .o\n");

    return 0;

}
