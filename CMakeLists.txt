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

# Windows
if(WIN32)

# Command-line version

set(PRAGUESKYMODEL_SRC
    thirdparty/miniz/miniz.c
    src/PragueSkyModel.cpp
    src/PragueSkyModelTestCli.cpp
)

add_executable(PragueSkyModelCli ${PRAGUESKYMODEL_SRC})

target_include_directories(PragueSkyModelCli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)
target_include_directories(PragueSkyModelCli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/miniz)

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

add_executable(PragueSkyModelGui ${PRAGUESKYMODELGUI_SRC})

target_include_directories(PragueSkyModelGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)
target_include_directories(PragueSkyModelGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/miniz)
target_include_directories(PragueSkyModelGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui)
target_include_directories(PragueSkyModelGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/backends)

target_link_libraries(PragueSkyModelGui PRIVATE d3d11.lib)

# Linux
else()

set(CMAKE_CXX_FLAGS "-DNDEBUG")

# Command-line version

find_package(TBB REQUIRED)

set(PRAGUESKYMODEL_SRC
    thirdparty/miniz/miniz.c
    src/PragueSkyModel.cpp
    src/PragueSkyModelTestCli.cpp
)

add_executable(PragueSkyModelCli ${PRAGUESKYMODEL_SRC})

target_include_directories(PragueSkyModelCli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)
target_include_directories(PragueSkyModelCli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/miniz)

target_link_libraries(PragueSkyModelCli PRIVATE TBB::tbb)

# GUI version

set(OpenGL_GL_PREFERENCE GLVND)

find_package(OpenGL REQUIRED)
find_package(SDL2 REQUIRED)

set(PRAGUESKYMODELGUI_SRC
    thirdparty/miniz/miniz.c
    thirdparty/imgui/imgui.cpp
    thirdparty/imgui/imgui_draw.cpp
    thirdparty/imgui/imgui_tables.cpp
    thirdparty/imgui/imgui_widgets.cpp
    thirdparty/imgui/backends/imgui_impl_sdl.cpp
    thirdparty/imgui/backends/imgui_impl_opengl2.cpp
    src/PragueSkyModel.cpp
    src/PragueSkyModelTestGui.cpp
)

add_executable(PragueSkyModelGui ${PRAGUESKYMODELGUI_SRC})

target_include_directories(PragueSkyModelGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)
target_include_directories(PragueSkyModelGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/miniz)
target_include_directories(PragueSkyModelGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui)
target_include_directories(PragueSkyModelGui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/imgui/backends)
target_include_directories(PragueSkyModelGui PRIVATE ${OPENGL_INCLUDE_DIRS})
target_include_directories(PragueSkyModelGui PRIVATE ${SDL2_INCLUDE_DIRS})

target_link_libraries(PragueSkyModelGui PRIVATE ${OPENGL_LIBRARIES})
target_link_libraries(PragueSkyModelGui PRIVATE ${SDL2_LIBRARIES})
target_link_libraries(PragueSkyModelGui PRIVATE TBB::tbb)
target_link_libraries(PragueSkyModelGui PRIVATE ${CMAKE_DL_LIBS})

endif()
