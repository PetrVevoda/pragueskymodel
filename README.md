# skymodel_intel

C++ implementation of [A Fitted Radiance and Attenuation Model for Realistic Atmospheres](https://cgg.mff.cuni.cz/publications/skymodel-2021/).

## Contents

- `src/PragueSkyModel.h`, `src/PragueSkyModel.cpp`
    - the implementation of the model
    - these two files are all you need to use the model in your code
- `src/PragueSkyModelTest.h`
    - a very simple example implementation of a render function that creates upfacing fisheye images of the sky
    - demonstrates how to use the model
- `src/PragueSkyModelTestCli.cpp`
    - a command-line user interface for `src/PragueSkyModelTest.h`
    - multiplatform
- `src/PragueSkyModelTestGui.cpp`
    - a graphical user interface for `src/PragueSkyModelTest.h`
    - currently works on Windows with DirectX 11 and on Linux with SDL 2 and OpenGL 2
- `thirdparty/tinyexr/tinyexr.h`, `thirdparty/miniz/miniz.h`, `thirdparty/miniz/miniz.c`
    - [Tiny OpenEXR image library](https://github.com/syoyo/tinyexr)
    - used by both user interfaces for saving results into EXR files
- `thirdparty/imgui`
    - [Dear ImGui library](https://github.com/ocornut/imgui)
    - used to build the graphical user interface
- `CMakeLists.txt`
    - cmake file for generating project files
- `README.md`
    - this readme
    
## Requirements

GUI version requires DirectX 11 on Windows and SDL 2 with OpenGL 2 on Linux.

## Usage

1. Clone the repository, run cmake and compile
    - tested on Windows 10 with Visual Studio 2022 and MSVC 19.31.31104.0 compiler
    - tested on Ubuntu 20.04 with GNU 9.3.0 compiler, SDL 2.0.10 and OpenGL 2.1
2. Download a dataset
    - the model is flexible and can work with various subsets of the original dataset
    - currently available:
        - [Full version (2.2 GB)](https://drive.google.com/file/d/19K96jEQmmqCeg8yjgZxj2awQj62lI50p/view?usp=sharing)
            - contains polarisation and entire range of ground albedos, observer altitudes, solar elevations, and visibilities as presented in the paper
        - [Ground-level version (102 MB)](https://drive.google.com/file/d/1Gk6OSHGpFx8HM3drHWykb3lDrtZXO4h7/view?usp=sharing)
            - contains all ground albedos, solar elevations, and visibilities, but only a single (zero) observer altitude, does not contain polarisation
3. Run the CLI version `PragueSkyModelCli.exe -dat <path_to_the_dataset>`
    - this will render a default configuration
    - use option `-h` or `--help` to display this list of available options:
        - `-alb` ... ground albedo, valid range [0, 1], default 0.5
        - `-alt` ... observer altitude, valid range [0, 15000] meters, default 0 meters
        - `-azi` ... solar azimuth, valid range [0, 360] degrees, default 0 degrees
        - `-cam` ... rendered view, 0 for up-facing fisheye, 1 for side-facing fisheye
        - `-dat` ... path to the dataset, default ".\PragueSkyModelDataset.dat"
        - `-ele` ... solar elevation, valid range [-4.2, 90] degrees, default 0 degrees
        - `-mod` ... what quantity to output, use 0 for sky radiance, 1 for sun radiance, 2 for polarisation and 3 for transmittance, default 0
        - `-out` ... output file, default "test.exr"
        - `-res` ... square output image resolution, default 512 pixels
        - `-vis` ... ground level visibility, valid range [20, 131.8] kilometers, default 59.4 kilometers
4. Or run the GUI version `PragueSkyModelGui.exe`
