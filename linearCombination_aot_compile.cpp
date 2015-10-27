// On linux, you can compile and run like so:
// g++ linearCombination_aot_compile.cpp -g -std=c++11 -I ./include -L ./bin -lHalide -lpthread -ldl -o lincombo_aot_generate
// LD_LIBRARY_PATH=./bin ./lincombo_aot_generate
// g++ linearCombination_aot_run.cpp lincombo_aot.o -lpthread -o lincombo_aot_run
// ./lincombo_aot_run

// On os x:
// g++ linearCombination_aot_compile.cpp -g -std=c++11 -I ./include -L ./bin -lHalide -o lincombo_aot_generate
// DYLD_LIBRARY_PATH=./bin ./lincombo_aot_generate
// g++ linearCombination_aot_run.cpp -g lincombo_aot.o -I ./include -I DarwinX86/pex_policy/10.1+1/include/ -I DarwinX86/daf_persistence/10.1+1/include/ -I DarwinX86/utils/10.1+1/include/ -I DarwinX86/daf_base/10.1+2/include/ -I DarwinX86/base/10.1+1/include/ -I DarwinX86/ndarray/10.1+2/include/ -I DarwinX86/pex_exceptions/10.1+1/include/ -I DarwinX86/eigen/3.2.0/include/ -I DarwinX86/afw/10.1+1/include -L ./bin -L DarwinX86/afw/10.1+1/lib -L DarwinX86/daf_base/10.1+2/lib/ -L DarwinX86/daf_persistence/10.1+1/lib/ -L DarwinX86/boost/1.55.0.1.lsst2+3/lib/ -lHalide -lafw -ldaf_base -ldaf_persistence -lboost_system `libpng-config --cflags --ldflags` -o lincombo_aot_run -std=c++11
// DYLD_LIBRARY_PATH=./bin:DarwinX86/afw/10.1+1/lib/:DarwinX86/daf_persistence/10.1+1/lib/:DarwinX86/daf_base/10.1+2/lib/:DarwinX86/boost/1.55.0.1.lsst2+3/lib/:DarwinX86/xpa/2.1.15.lsst2/lib/:DarwinX86/pex_policy/10.1+1/lib/:DarwinX86/pex_logging/10.1+1/lib/:DarwinX86/utils/10.1+1/lib/:DarwinX86/pex_exceptions/10.1+1/lib/:DarwinX86/base/10.1+1/lib/ ./lincombo_aot_run


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

    ImageParam polynomialCoefficients(type_of<float>(), 2); //array with dimension [5][10]
    //Kernel parameters array with dimensions [5][3]
    //second dimension: sigmaX, sigmaY, theta
    ImageParam kerParams(type_of<float>(), 2);

    //Kernel has dimensions (boundingBox*2 + 1) x (boundingBox*2 + 1)
    int boundingBox = 2; 
    Var x, y, y0, yi;
    float pi = 3.14159265359f;

    Func polynomials[num_kernels];
    for(int k = 0; k < num_kernels; k++){
        polynomials[k](x, y) = polynomialCoefficients(k,0) + 
            polynomialCoefficients(k,1)*x + polynomialCoefficients(k,2)*y +
            polynomialCoefficients(k,3)*x*x + polynomialCoefficients(k,4)*x*y + 
            polynomialCoefficients(k,5)*y*y + polynomialCoefficients(k,6)*x*x*x +
            polynomialCoefficients(k,7)*x*x*y + polynomialCoefficients(k,8)*x*y*y
            + polynomialCoefficients(k,9)*y*y*y;
    }

    Func kernels[num_kernels];
    Var i,j;
    for(int k = 0; k < num_kernels; k++){
        kernels[k](i, j) = exp(-((i*cos(kerParams(k,2)) + j*sin(kerParams(k,2)))*
                    (i*cos(kerParams(k,2)) + j*sin(kerParams(k,2))))
                    / (2*kerParams(k,0)*kerParams(k,0))
                    -((j*cos(kerParams(k,2)) - i*sin(kerParams(k,2)))*
                    (j*cos(kerParams(k,2)) - i*sin(kerParams(k,2))))
                    / (2*kerParams(k,1)*kerParams(k,1))) /
                    (2.0f*pi*kerParams(k,0)*kerParams(k,1));
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

    Expr curKernelVal = 0.0f;
    for(int i = -boundingBox; i <= boundingBox; i++){
        for(int j = -boundingBox; j <= boundingBox; j++){
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

    //Evaluate image, mask, and variance planes concurrently using a tuple
    Func combined_output ("combined_output");
    combined_output(x, y) = Tuple(blur_image_help, blur_variance_help, maskOutHelp);

    // Split the y coordinate of the output into strips of 32 scanlines:
    combined_output.split(y, y0, yi, 32);
    // Compute the strips using a thread pool and a task queue.
    combined_output.parallel(y0);
    // Vectorize across x by a factor of eight.
    combined_output.vectorize(x, 8);

    std::vector<Argument> args = {image, variance, mask, polynomialCoefficients, kerParams};
    combined_output.compile_to_file("lincombo_aot", args);

    printf("Halide pipeline compiled, but not yet run.\n");

    return 0;

}


//    //Five 3rd degree polynomials which will be used as spatially varying
//    //coefficients in the linear combination of the five gaussian basis kernels
//    Func polynomial1 ("polynomial1");
//    polynomial1(x, y) = 0.1f + 0.002f*x + 0.003f*y + 0.4f*x*x + 0.5f*x*y
//                     + 0.6f*y*y +  0.0007f*x*x*x + 0.0008f*x*x*y + 0.0009f*x*y*y
//                     + 0.00011f*y*y*y;
//
//    Func polynomial2 ("polynomial2");
//    polynomial2(x, y) = 1.1f + 1.002f*x + 1.003f*y + 1.4f*x*x + 1.5f*x*y
//                     + 1.6f*y*y +  1.0007f*x*x*x + 1.0008f*x*x*y + 1.0009f*x*y*y
//                     + 1.00011f*y*y*y;
//
//    Func polynomial3 ("polynomial3");
//    polynomial3(x, y) = 2.1f + 2.002f*x + 2.003f*y + 2.4f*x*x + 2.5f*x*y
//                     + 2.6f*y*y +  2.0007f*x*x*x + 2.0008f*x*x*y + 2.0009f*x*y*y
//                     + 2.00011f*y*y*y;
//
//    Func polynomial4 ("polynomial4");
//    polynomial4(x, y) = 3.1f + 3.002f*x + 3.003f*y + 3.4f*x*x + 3.5f*x*y
//                     + 3.6f*y*y +  3.0007f*x*x*x + 3.0008f*x*x*y + 3.0009f*x*y*y
//                     + 3.00011f*y*y*y;
//
//    Func polynomial5 ("polynomial5");
//    polynomial5(x, y) = 4.1f + 4.002f*x + 4.003f*y + 4.4f*x*x + 4.5f*x*y
//                     + 4.6f*y*y +  4.0007f*x*x*x + 4.0008f*x*x*y + 4.0009f*x*y*y
//                     + 4.00011f*y*y*y;


    //5 Guassian basis kernels
    //Kernel #1
//    Func kernel1;
//    float sigmaX1 = 2.0f;
//    float sigmaY1 = 2.0f;
//    float theta1 = 0.0f; //rotation of sigmaX axis
//    kernel1(x, y) = exp(-((x*cos(theta1) + y*sin(theta1))*(x*cos(theta1) + y*sin(theta1)))
//                    /(2*sigmaX1*sigmaX1)
//                    -((y*cos(theta1) - x*sin(theta1))*(y*cos(theta1) - x*sin(theta1)))
//                    /(2*sigmaY1*sigmaY1)) / (2.0f*pi*sigmaX1*sigmaY1);
//
//    //Kernel #2
//    Func kernel2;
//    float sigmaX2 = 0.5f;
//    float sigmaY2 = 4.0f;
//    float theta2 = 0.0f; //rotation of sigmaX axis
//    kernel2(x, y) = exp(-((x*cos(theta2) + y*sin(theta2))*(x*cos(theta2) + y*sin(theta2)))
//                    /(2*sigmaX2*sigmaX2)
//                    -((y*cos(theta2) - x*sin(theta2))*(y*cos(theta2) - x*sin(theta2)))
//                    /(2*sigmaY2*sigmaY2)) / (2.0f*pi*sigmaX2*sigmaY2);
//
//    //Kernel #3
//    Func kernel3;
//    float sigmaX3 = 0.5f;
//    float sigmaY3 = 4.0f;
//    float theta3 = 3.14159f/4; //rotation of sigmaX axis
//    kernel3(x, y) = exp(-((x*cos(theta3) + y*sin(theta3))*(x*cos(theta3) + y*sin(theta3)))
//                    /(2*sigmaX3*sigmaX3)
//                    -((y*cos(theta3) - x*sin(theta3))*(y*cos(theta3) - x*sin(theta3)))
//                    /(2*sigmaY3*sigmaY3)) / (2.0f*pi*sigmaX3*sigmaY3);
//    //Kernel #4
//    Func kernel4;
//    float sigmaX4 = 0.5f;
//    float sigmaY4 = 4.0f;
//    float theta4 = 3.14159f/2; //rotation of sigmaX axis
//    kernel4(x, y) = exp(-((x*cos(theta4) + y*sin(theta4))*(x*cos(theta4) + y*sin(theta4)))
//                    /(2*sigmaX4*sigmaX4)
//                    -((y*cos(theta4) - x*sin(theta4))*(y*cos(theta4) - x*sin(theta4)))
//                    /(2*sigmaY4*sigmaY4)) / (2.0f*pi*sigmaX4*sigmaY4);
//
//
//    //Kernel #5
//    Func kernel5;
//    float sigmaX5 = 4.0f;
//    float sigmaY5 = 4.0f;
//    float theta5 = 0.0; //rotation of sigmaX axis
//    kernel5(x, y) = (exp(-((x*cos(theta5) +y*sin(theta5))*(x*cos(theta5) +y*sin(theta5)))
//                    /(2*sigmaX5*sigmaX5)))
//                    *(exp(-((y*cos(theta5) - x*sin(theta5))*(y*cos(theta5) - x*sin(theta5)))
//                    /(2*sigmaY5*sigmaY5)) / (2.0f*pi*sigmaX5*sigmaY5));
