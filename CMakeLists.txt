## Copyright 2022 Intel Corporation
## SPDX-License-Identifier: Apache-2.0

## Global CMake options ##

cmake_minimum_required(VERSION 3.1)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(PragueSkyModel)

# Command-line version

set(PRAGUESKYMODEL_SRC
    thirdparty/miniz/miniz.c
    src/PragueSkyModel.cpp
    src/PragueSkyModelTestCli.cpp
)

add_executable(PragueSkyModelDemoCli ${PRAGUESKYMODEL_SRC})

target_include_directories(PragueSkyModelDemoCli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)
target_include_directories(PragueSkyModelDemoCli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/miniz)

# GUI version

find_library(D3D11_LIBRARY d3d11.lib)

set(PRAGUESKYMODELGUI_SRC
    thirdparty/miniz/miniz.c
    thirdparty/imgui/imgui.cpp
    thirdparty/imgui/imgui_draw.cpp
    thirdparty/imgui/imgui_tables.cpp
    thirdparty/imgui/imgui_widgets.cpp
    thirdparty/imgui/backends/imgui_impl_dx11.cpp
    thirdparty/imgui/backends/imgui_impl_win32.cpp
    src/PragueSkyModel.cpp
    src/PragueSkyModelTestGui.cpp
)

add_executable(PragueSkyModelDemoGui ${PRAGUESKYMODELGUI_SRC})

target_include_directories(PragueSkyModelDemoGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)
target_include_directories(PragueSkyModelDemoGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/miniz)
target_include_directories(PragueSkyModelDemoGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui)
target_include_directories(PragueSkyModelDemoGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/backends)

target_link_libraries(PragueSkyModelDemoGui PRIVATE d3d11.lib)