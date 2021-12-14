#include <array>
#include <iostream>
#include <limits>
#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>

// For saving the result to an EXR file.
#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"

#include "PragueSkyModel.h"


using Vector3  = PragueSkyModel::Vector3;

// We use 11-channel spectrum for querying the model. The wavelengths samples are placed at centers of
// spectral bins used by the model.
constexpr int    SPECTRUM_CHANNELS      = 11;
constexpr double SPECTRUM_STEP          = 40;
using Spectrum                          = std::array<double, SPECTRUM_CHANNELS>;
constexpr Spectrum SPECTRUM_WAVELENGTHS = { 340.0, 380.0, 420.0, 460.0, 500.0, 540.0,
                                            580.0, 620.0, 660.0, 700.0, 740.0 };


/////////////////////////////////////////////////////////////////////////////////////
// Command line
/////////////////////////////////////////////////////////////////////////////////////

/// Checks if given command-line option was used.
bool cmdOptionExists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}

/// Tries to get string value of given command-line option. Returns default value if not successful.
std::string getStringCmdOption(char**             begin,
                               char**             end,
                               const std::string& option,
                               const std::string& defaultValue) {
    char** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
        return std::string(*itr);
    } else {
        return defaultValue;
    }
}

/// Tries to get double value of given command-line option. Returns default value if not successful.
double getDoubleCmdOption(char** begin, char** end, const std::string& option, const double defaultValue) {
    char** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
        try {
            return std::stod(*itr);
        } catch (std::invalid_argument const& ex) {
            std::cout << "Warning: Invalid value of " << option << " argument, using default " << defaultValue
                      << "\n";
        }
    }
    return defaultValue;
}

/// Tries to get integer value of given command-line option. Returns default value if not successful.
int getIntCmdOption(char** begin, char** end, const std::string& option, const int defaultValue) {
    char** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
        try {
            return std::stoi(*itr);
        } catch (std::invalid_argument const& ex) {
            std::cout << "Warning: Invalid value of " << option << " argument, using default " << defaultValue
                      << "\n";
        }
    }
    return defaultValue;
}


/////////////////////////////////////////////////////////////////////////////////////
// Conversion functions
/////////////////////////////////////////////////////////////////////////////////////

double degreesToRadians(const double degrees) {
    return degrees * M_PI / 180.0;
}

/// Computes direction corresponding to given pixel coordinates in fisheye projection.
Vector3 pixelToDirection(int x, int y, int resolution) {
    // Make circular image area in center of image.
    const double radius  = resolution / 2;
    const double scaledx = (x - radius) / radius;
    const double scaledy = (y - radius) / radius;
    const double denom   = scaledx * scaledx + scaledy * scaledy + 1;

    if (denom > 2.0) {
        // Outside image area.
        return Vector3();
    } else {
        // Stereographic mapping.
        return Vector3(2.0 * scaledx / denom, 2.0 * scaledy / denom, -(denom - 2.0) / denom);
    }
}

// Spectral response table used for converting spectrum to XYZ.
const Vector3 SPECTRAL_RESPONSE[] = {
    Vector3(0.000129900000f, 0.000003917000f, 0.000606100000f),
    Vector3(0.000232100000f, 0.000006965000f, 0.001086000000f),
    Vector3(0.000414900000f, 0.000012390000f, 0.001946000000f),
    Vector3(0.000741600000f, 0.000022020000f, 0.003486000000f),
    Vector3(0.001368000000f, 0.000039000000f, 0.006450001000f),
    Vector3(0.002236000000f, 0.000064000000f, 0.010549990000f),
    Vector3(0.004243000000f, 0.000120000000f, 0.020050010000f),
    Vector3(0.007650000000f, 0.000217000000f, 0.036210000000f),
    Vector3(0.014310000000f, 0.000396000000f, 0.067850010000f),
    Vector3(0.023190000000f, 0.000640000000f, 0.110200000000f),
    Vector3(0.043510000000f, 0.001210000000f, 0.207400000000f),
    Vector3(0.077630000000f, 0.002180000000f, 0.371300000000f),
    Vector3(0.134380000000f, 0.004000000000f, 0.645600000000f),
    Vector3(0.214770000000f, 0.007300000000f, 1.039050100000f),
    Vector3(0.283900000000f, 0.011600000000f, 1.385600000000f),
    Vector3(0.328500000000f, 0.016840000000f, 1.622960000000f),
    Vector3(0.348280000000f, 0.023000000000f, 1.747060000000f),
    Vector3(0.348060000000f, 0.029800000000f, 1.782600000000f),
    Vector3(0.336200000000f, 0.038000000000f, 1.772110000000f),
    Vector3(0.318700000000f, 0.048000000000f, 1.744100000000f),
    Vector3(0.290800000000f, 0.060000000000f, 1.669200000000f),
    Vector3(0.251100000000f, 0.073900000000f, 1.528100000000f),
    Vector3(0.195360000000f, 0.090980000000f, 1.287640000000f),
    Vector3(0.142100000000f, 0.112600000000f, 1.041900000000f),
    Vector3(0.095640000000f, 0.139020000000f, 0.812950100000f),
    Vector3(0.057950010000f, 0.169300000000f, 0.616200000000f),
    Vector3(0.032010000000f, 0.208020000000f, 0.465180000000f),
    Vector3(0.014700000000f, 0.258600000000f, 0.353300000000f),
    Vector3(0.004900000000f, 0.323000000000f, 0.272000000000f),
    Vector3(0.002400000000f, 0.407300000000f, 0.212300000000f),
    Vector3(0.009300000000f, 0.503000000000f, 0.158200000000f),
    Vector3(0.029100000000f, 0.608200000000f, 0.111700000000f),
    Vector3(0.063270000000f, 0.710000000000f, 0.078249990000f),
    Vector3(0.109600000000f, 0.793200000000f, 0.057250010000f),
    Vector3(0.165500000000f, 0.862000000000f, 0.042160000000f),
    Vector3(0.225749900000f, 0.914850100000f, 0.029840000000f),
    Vector3(0.290400000000f, 0.954000000000f, 0.020300000000f),
    Vector3(0.359700000000f, 0.980300000000f, 0.013400000000f),
    Vector3(0.433449900000f, 0.994950100000f, 0.008749999000f),
    Vector3(0.512050100000f, 1.000000000000f, 0.005749999000f),
    Vector3(0.594500000000f, 0.995000000000f, 0.003900000000f),
    Vector3(0.678400000000f, 0.978600000000f, 0.002749999000f),
    Vector3(0.762100000000f, 0.952000000000f, 0.002100000000f),
    Vector3(0.842500000000f, 0.915400000000f, 0.001800000000f),
    Vector3(0.916300000000f, 0.870000000000f, 0.001650001000f),
    Vector3(0.978600000000f, 0.816300000000f, 0.001400000000f),
    Vector3(1.026300000000f, 0.757000000000f, 0.001100000000f),
    Vector3(1.056700000000f, 0.694900000000f, 0.001000000000f),
    Vector3(1.062200000000f, 0.631000000000f, 0.000800000000f),
    Vector3(1.045600000000f, 0.566800000000f, 0.000600000000f),
    Vector3(1.002600000000f, 0.503000000000f, 0.000340000000f),
    Vector3(0.938400000000f, 0.441200000000f, 0.000240000000f),
    Vector3(0.854449900000f, 0.381000000000f, 0.000190000000f),
    Vector3(0.751400000000f, 0.321000000000f, 0.000100000000f),
    Vector3(0.642400000000f, 0.265000000000f, 0.000049999990f),
    Vector3(0.541900000000f, 0.217000000000f, 0.000030000000f),
    Vector3(0.447900000000f, 0.175000000000f, 0.000020000000f),
    Vector3(0.360800000000f, 0.138200000000f, 0.000010000000f),
    Vector3(0.283500000000f, 0.107000000000f, 0.000000000000f),
    Vector3(0.218700000000f, 0.081600000000f, 0.000000000000f),
    Vector3(0.164900000000f, 0.061000000000f, 0.000000000000f),
    Vector3(0.121200000000f, 0.044580000000f, 0.000000000000f),
    Vector3(0.087400000000f, 0.032000000000f, 0.000000000000f),
    Vector3(0.063600000000f, 0.023200000000f, 0.000000000000f),
    Vector3(0.046770000000f, 0.017000000000f, 0.000000000000f),
    Vector3(0.032900000000f, 0.011920000000f, 0.000000000000f),
    Vector3(0.022700000000f, 0.008210000000f, 0.000000000000f),
    Vector3(0.015840000000f, 0.005723000000f, 0.000000000000f),
    Vector3(0.011359160000f, 0.004102000000f, 0.000000000000f),
    Vector3(0.008110916000f, 0.002929000000f, 0.000000000000f),
    Vector3(0.005790346000f, 0.002091000000f, 0.000000000000f),
    Vector3(0.004109457000f, 0.001484000000f, 0.000000000000f),
    Vector3(0.002899327000f, 0.001047000000f, 0.000000000000f),
    Vector3(0.002049190000f, 0.000740000000f, 0.000000000000f),
    Vector3(0.001439971000f, 0.000520000000f, 0.000000000000f),
    Vector3(0.000999949300f, 0.000361100000f, 0.000000000000f),
    Vector3(0.000690078600f, 0.000249200000f, 0.000000000000f),
    Vector3(0.000476021300f, 0.000171900000f, 0.000000000000f),
    Vector3(0.000332301100f, 0.000120000000f, 0.000000000000f),
    Vector3(0.000234826100f, 0.000084800000f, 0.000000000000f),
    Vector3(0.000166150500f, 0.000060000000f, 0.000000000000f),
    Vector3(0.000117413000f, 0.000042400000f, 0.000000000000f),
    Vector3(0.000083075270f, 0.000030000000f, 0.000000000000f),
    Vector3(0.000058706520f, 0.000021200000f, 0.000000000000f),
    Vector3(0.000041509940f, 0.000014990000f, 0.000000000000f),
    Vector3(0.000029353260f, 0.000010600000f, 0.000000000000f),
    Vector3(0.000020673830f, 0.000007465700f, 0.000000000000f),
    Vector3(0.000014559770f, 0.000005257800f, 0.000000000000f),
    Vector3(0.000010253980f, 0.000003702900f, 0.000000000000f),
    Vector3(0.000007221456f, 0.000002607800f, 0.000000000000f),
    Vector3(0.000005085868f, 0.000001836600f, 0.000000000000f),
    Vector3(0.000003581652f, 0.000001293400f, 0.000000000000f),
    Vector3(0.000002522525f, 0.000000910930f, 0.000000000000f),
    Vector3(0.000001776509f, 0.000000641530f, 0.000000000000f),
    Vector3(0.000001251141f, 0.000000451810f, 0.000000000000f),
};
constexpr float SPECTRAL_RESPONSE_START = 360.f;
constexpr float SPECTRAL_RESPONSE_STEP  = 5.f;

/// Converts given spectrum to sRGB.
Vector3 spectrumToRGB(const Spectrum& spectrum) {
    // Spectrum to XYZ
    Vector3 xyz = Vector3();
    for (int wl = 0; wl < SPECTRUM_CHANNELS; wl++) {
        const int responseIdx = (SPECTRUM_WAVELENGTHS[wl] - SPECTRAL_RESPONSE_START) / SPECTRAL_RESPONSE_STEP;
        xyz                   = xyz + SPECTRAL_RESPONSE[responseIdx] * spectrum[wl];
    }
    xyz = xyz * SPECTRUM_STEP;

    // XYZ to sRGB
    Vector3 rgb = Vector3();
    rgb.x       = 3.2404542 * xyz.x - 1.5371385 * xyz.y - 0.4985314 * xyz.z;
    rgb.y       = -0.9692660 * xyz.x + 1.8760108 * xyz.y + 0.0415560 * xyz.z;
    rgb.z       = 0.0556434 * xyz.x - 0.2040259 * xyz.y + 1.0572252 * xyz.z;

    return rgb;
}



/////////////////////////////////////////////////////////////////////////////////////
// EXR
/////////////////////////////////////////////////////////////////////////////////////

// Function for saving image data into an EXR file. Code from TinyEXR.
void saveEXR(const std::vector<float>& rgb, const int width, const int height, const std::string outfilename) {
    // Header

    EXRHeader header;
    InitEXRHeader(&header);

    header.num_channels = 3;
    header.channels     = (EXRChannelInfo*)malloc(sizeof(EXRChannelInfo) * header.num_channels);
    // Must be (A)BGR order, since most of EXR viewers expect this channel order.
    strncpy(header.channels[0].name, "B", 255);
    header.channels[0].name[strlen("B")] = '\0';
    strncpy(header.channels[1].name, "G", 255);
    header.channels[1].name[strlen("G")] = '\0';
    strncpy(header.channels[2].name, "R", 255);
    header.channels[2].name[strlen("R")] = '\0';

    header.pixel_types           = (int*)malloc(sizeof(int) * header.num_channels);
    header.requested_pixel_types = (int*)malloc(sizeof(int) * header.num_channels);
    for (int i = 0; i < header.num_channels; i++) {
        header.pixel_types[i]           = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
        header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of output image to be stored
                                                                   // in .EXR
    }

    // Image

    EXRImage image;
    InitEXRImage(&image);

    image.num_channels = 3;

    const size_t       pixelCount = size_t(width) * height;
    std::vector<float> images[3];
    images[0].resize(pixelCount);
    images[1].resize(pixelCount);
    images[2].resize(pixelCount);

    // Split RGBRGBRGB... into R, G and B layer
    for (size_t i = 0; i < pixelCount; i++) {
        images[0][i] = rgb[3 * i + 0];
        images[1][i] = rgb[3 * i + 1];
        images[2][i] = rgb[3 * i + 2];
    }

    float* imagePtr[3];
    imagePtr[0] = &(images[2].at(0)); // B
    imagePtr[1] = &(images[1].at(0)); // G
    imagePtr[2] = &(images[0].at(0)); // R

    image.images = (unsigned char**)imagePtr;
    image.width  = width;
    image.height = height;

    // Write

    const char* err = nullptr;
    int         ret = SaveEXRImageToFile(&image, &header, outfilename.c_str(), &err);
    if (ret != TINYEXR_SUCCESS) {
        const std::string message(std::string("Saving EXR failed - ") + std::string(err));
        FreeEXRErrorMessage(err); // Frees buffer for an error message
        throw std::exception(message.c_str());
    } else {
        free(header.channels);
        free(header.pixel_types);
        free(header.requested_pixel_types);
    }
}



/////////////////////////////////////////////////////////////////////////////////////
// Main
/////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
    // Print help if user asked for it and quit.
    if (cmdOptionExists(argv, argv + argc, "-h") || cmdOptionExists(argv, argv + argc, "--help")) {
        std::cout << "Usage: " << argv[0] << "\n";
        std::cout << "Optional arguments:\n";
        std::cout << "   -alb ... ground albedo, valid range [0, 1], default 0.5\n";
        std::cout << "   -alt ... observer altitude, valid range [0, 15000] meters, default 0 meters\n";
        std::cout << "   -azi ... solar azimuth, valid range [0, 360] degrees, default 0 degrees\n";
        std::cout << "   -dat ... path to the dataset, default \".\\PragueSkyModelDataset.dat\"\n";
        std::cout << "   -ele ... solar elevation, valid range [-4.2, 90] degrees, default 0 degrees\n";
        std::cout << "   -mod ... what quantity to output, use 0 for sky radiance, 1 for sun radiance, 2 for "
                     "polarisation and 3 for transmittance, default 0\n";
        std::cout << "   -out ... output file, default \"test.exr\"\n";
        std::cout << "   -res ... square output image resolution, default 512 pixels\n";
        std::cout << "   -vis ... ground level visibility, valid range [20, 131.8] kilometers, default 59.4 "
                     "kilometers\n";
        return 1;
    }

    // Read all command-line options.
    const double      albedo     = getDoubleCmdOption(argv, argv + argc, "-alb", 0.5);
    const double      altitude   = getDoubleCmdOption(argv, argv + argc, "-alt", 0.0);
    const double      azimuth    = degreesToRadians(getDoubleCmdOption(argv, argv + argc, "-azi", 0.0));
    const std::string dataset    = getStringCmdOption(argv, argv + argc, "-dat", "PragueSkyModelDataset.dat");
    const double      elevation  = degreesToRadians(getDoubleCmdOption(argv, argv + argc, "-ele", 0.0));
    const int         mode       = getIntCmdOption(argv, argv + argc, "-mod", 0);
    const std::string outputFile = getStringCmdOption(argv, argv + argc, "-out", "test.exr");
    const int         resolution = getIntCmdOption(argv, argv + argc, "-res", 512);
    const double      visibility = getDoubleCmdOption(argv, argv + argc, "-vis", 59.4);
    
    // We are viewing the sky from 'altitude' meters above the origin.
    const Vector3 viewPoint = Vector3(0.0, 0.0, altitude);

    // Buffer for saving the resulting image to an EXR file.
    std::vector<float> result;
    result.resize(size_t(resolution) * resolution * 3);

    // The model can throw exceptions, therefore try.
    try {
        // Initialize the model with the given dataset file.
        PragueSkyModel skyModel = PragueSkyModel(dataset);

        for (int x = 0; x < resolution; x++) {
            for (int y = 0; y < resolution; y++) {
                // For each pixel of the rendered image get the corresponding direction in fisheye projection.
                const Vector3 viewDir = pixelToDirection(x, y, resolution);

                // If the pixel lies outside the upper hemisphere, the direction will be zero. Such a pixel is
                // painted black.
                if (viewDir.isZero()) {
                    result[(size_t(x) * resolution + y) * 3]     = 0.0;
                    result[(size_t(x) * resolution + y) * 3 + 1] = 0.0;
                    result[(size_t(x) * resolution + y) * 3 + 2] = 0.0;
                    continue;
                }

                // Get internal model parameters for the desired configuration.
                const PragueSkyModel::Parameters params =
                    skyModel.computeParameters(viewPoint, viewDir, elevation, azimuth, visibility, albedo);

                // Based on the selected mode compute spectral sky radiance, sun radiance, polarisation or
                // transmittance.
                Spectrum spectrum;
                for (int wl = 0; wl < SPECTRUM_CHANNELS; wl++) {
                    switch (mode) {
                    case 1:
                        spectrum[wl] = skyModel.sunRadiance(params, SPECTRUM_WAVELENGTHS[wl]);
                        break;
                    case 2:
                        spectrum[wl] = std::abs(skyModel.polarisation(params, SPECTRUM_WAVELENGTHS[wl]));
                        break;
                    case 3:
                        spectrum[wl] = skyModel.transmittance(params,
                                                              SPECTRUM_WAVELENGTHS[wl],
                                                              std::numeric_limits<double>::max());
                        break;
                    default:
                        spectrum[wl] = skyModel.skyRadiance(params, SPECTRUM_WAVELENGTHS[wl]);
                        break;
                    }
                }

                // Convert the spectral quantity to sRGB and store it in the result buffer.
                const Vector3 rgb = spectrumToRGB(spectrum);
                result[(size_t(x) * resolution + y) * 3] = rgb.x;
                result[(size_t(x) * resolution + y) * 3 + 1] = rgb.y;
                result[(size_t(x) * resolution + y) * 3 + 2] = rgb.z;
            }
        }

        // Save the result buffer into an EXR file.
        saveEXR(result, resolution, resolution, outputFile);

		std::cout << "Done\n";
		return 0;
	}
	catch (std::exception& e) {
		std::cout << "Error: " << e.what() << "\n";
        return 1;
	}
}