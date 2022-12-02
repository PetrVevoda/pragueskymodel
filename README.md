# Prague Sky Model

This a C++ implementation of a sky model developed at Charles University in Prague and described in

> Wilkie et al. [A fitted radiance and attenuation model for realistic atmospheres](https://cgg.mff.cuni.cz/publications/skymodel-2021/). ACM Transactions on Graphics (Proceedings of SIGGRAPH 2021), 40(4), 2021.

and in

> Vévoda et al. [A Wide Spectral Range Sky Radiance Model](https://cgg.mff.cuni.cz/publications/infrared-skymodel-2022/). Computer Graphics Forum, 41(7), 2022.

The model is provided together with a simple renderer that demonstrates its usage and allows direct generation of sky images for various parameters.

## Contents

- `src/PragueSkyModel.h`, `src/PragueSkyModel.cpp`
    - the implementation of the model
    - these two files are all you need to use the model in your code
- `src/PragueSkyModelTest.h`
    - a very simple example implementation of a render function that creates fisheye images of the sky
    - demonstrates how to use the model
- `src/PragueSkyModelTestCli.cpp`
    - a command-line user interface for `src/PragueSkyModelTest.h`
    - multiplatform
- `src/PragueSkyModelTestGui.cpp`
    - a graphical user interface for `src/PragueSkyModelTest.h`
    - currently works on Windows with DirectX 11 and on Linux with SDL 2 and OpenGL 2
- `thirdparty/tinyexr`, `thirdparty/miniz`
    - [Tiny OpenEXR image library](https://github.com/syoyo/tinyexr)
    - used by both user interfaces for saving results into EXR files
- `thirdparty/imgui`
    - [Dear ImGui library](https://github.com/ocornut/imgui) and [imgui-filebrowser](https://github.com/AirGuanZ/imgui-filebrowser)
    - used to build the graphical user interface
- `CMakeLists.txt`
    - cmake file for generating project files
- `LICENSE.txt`
    - Apache-2.0 license file
- `README.md`
    - this readme
    
## Requirements

C++ 17 is required. GUI version also requires DirectX 11 on Windows and SDL 2 with OpenGL 2 on Linux. Both CLI and GUI versions require TBB on Linux.

## Usage

1. Clone the repository, run cmake and compile
    - tested on Windows 10 with Visual Studio 2022 and MSVC 19.31.31104.0 compiler
    - tested on Ubuntu 20.04 with GNU 9.3.0 compiler, SDL 2.0.10 and OpenGL 2.1
2. Download a dataset
    - the model is flexible and can work with various versions of the dataset
    - currently available:
        - [Full version (2.2 GB)](https://drive.google.com/file/d/19K96jEQmmqCeg8yjgZxj2awQj62lI50p/view?usp=sharing)
            - a full version of the dataset as presented in the first paper, contains the entire range of visibilities, solar elevations, observer altitudes, and ground albedos, includes polarisation
        - [Ground-level version (103 MB)](https://drive.google.com/file/d/1IflyFZTJxC_N298yXq_2GK4ycIsVJZk6/view?usp=sharing)
            - a smaller version of the dataset, contains only a single (zero) observer altitude, does not include polarisation
        - [SWIR version (547 MB)](https://drive.google.com/file/d/1ZOizQCN6tH39JEwyX8KvAj7WEdX-EqJl/view?usp=sharing)
            - a wide spectral range version of the dataset as presented in the second paper, contains 55 (instead of 11) wavelength channels, but only a single (zero) observer altitude, includes polarisation

            |                         | Full               | Ground-level       | SWIR                |
            | ----------------------- |:------------------:|:------------------:|:-------------------:|
            | **visibilities**        | all                | all                | all                 |
            | **solar elevations**    | all                | all                | all                 |
            | **observer altitudes**  | all                | just one (0 m)     | just one (0 m)      |
            | **ground albedos**      | all                | all                | all                 |
            | **wavelength channels** | 11 (320 - 760 nm)  | 11 (320 - 760 nm)  | 55 (280 - 2480 nm)  |
            | **transmittance**       | yes                | yes                | yes                 |
            | **polarisation**        | yes                | no                 | yes                 |

3. Run the CLI version `PragueSkyModelCli -dat <path_to_the_dataset>`
    - this will render a default configuration
    - use option `-h` or `--help` to display this list of available options:
        - `-alb` ... ground albedo, valid range [0, 1], default 0.5
        - `-alt` ... observer altitude, valid range [0, 15000] meters, default 0 meters
        - `-azi` ... solar azimuth, valid range [0, 360] degrees, default 0 degrees
        - `-cam` ... rendered view, 0 for up-facing fisheye, 1 for side-facing fisheye
        - `-chn` ... output channel, 0 for sRGB with visible range, 1 - 55 for individual channels (from 280 to 2480 nm with 40 nm steps)
        - `-dat` ... path to the dataset, default "PragueSkyModelDataset.dat"
        - `-ele` ... solar elevation, valid range [-4.2, 90] degrees, default 0 degrees
        - `-mod` ... what quantity to output, use 0 for sky radiance, 1 for sun radiance, 2 for polarisation and 3 for transmittance, default 0
        - `-out` ... output file, default "test.exr"
        - `-res` ... square output image resolution, default 512 pixels
        - `-vis` ... ground level visibility, valid range [20, 131.8] kilometers, default 59.4 kilometers
4. Or run the GUI version using `PragueSkyModelGui`

## Acknowledgments

This implementation was supported by [Chaos Czech](https://corona-renderer.com/) and the Intel Graphics and Visualization Institute at Charles University.
