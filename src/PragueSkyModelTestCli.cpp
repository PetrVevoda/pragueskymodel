// For saving the result to an EXR file.
#pragma warning(push)
#pragma warning(disable : 26495)
#pragma warning(disable : 26451)
#pragma warning(disable : 6386)
#pragma warning(disable : 6387)
#pragma warning(disable : 6385)
#pragma warning(disable : 6011)
#pragma warning(disable : 6001)
#pragma warning(disable : 4018)
#define TINYEXR_IMPLEMENTATION
#include <tinyexr/tinyexr.h>
#pragma warning(pop)

// For CLI
#include <iostream>

#include "PragueSkyModelTest.h"


/////////////////////////////////////////////////////////////////////////////////////
// Command line helper functions
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
        } catch (std::invalid_argument) {
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
        } catch (std::invalid_argument) {
            std::cout << "Warning: Invalid value of " << option << " argument, using default " << defaultValue
                      << "\n";
        }
    }
    return defaultValue;
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
        std::cout << "   -cam ... rendered view, 0 for up-facing fisheye, 1 for side-facing fisheye\n";
        std::cout << "   -chn ... output channel, 0 for sRGB with visible range, 1 - 55 for individual channels (280:40:2480 nm)\n";
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
    const int         channel    = std::clamp(getIntCmdOption(argv, argv + argc, "-chn", 0), 0, SPECTRUM_CHANNELS);
    const std::string dataset    = getStringCmdOption(argv, argv + argc, "-dat", "PragueSkyModelDataset.dat");
    const double      elevation  = degreesToRadians(getDoubleCmdOption(argv, argv + argc, "-ele", 0.0));
    const Mode        mode       = Mode(getIntCmdOption(argv, argv + argc, "-mod", 0));
    const std::string outputFile = getStringCmdOption(argv, argv + argc, "-out", "test.exr");
    const int         resolution = getIntCmdOption(argv, argv + argc, "-res", 512);
    const View        view       = View(getIntCmdOption(argv, argv + argc, "-cam", 0));
    const double      visibility = getDoubleCmdOption(argv, argv + argc, "-vis", 59.4);

    std::vector<std::vector<float>> result;

    // The model can throw exceptions, therefore try.
    try {
        // Initialize the model with the given dataset file.
        PragueSkyModel skyModel;
        skyModel.initialize(dataset, visibility);

        // Render sky image according to the given configuration.
        render(skyModel, albedo, altitude, azimuth, elevation, mode, resolution, view, visibility, result);

        // Save the result buffer into an EXR file.
        const char* err = nullptr;
        int ret = SaveEXR(result[channel].data(), resolution, resolution, channel == 0 ? 3 : 1, 0, outputFile.c_str(), &err);
        if (ret != TINYEXR_SUCCESS) {
            const std::string message(std::string("Saving EXR failed - ") + std::string(err));
            FreeEXRErrorMessage(err); // Frees buffer for an error message
            throw message.c_str();
        }

        std::cout << "Done\n";
        return 0;
    } catch (std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }
}