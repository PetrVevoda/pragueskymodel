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

project(pragueskymodel)


set(PRAGUESKYMODEL_SRC
    miniz.c
    PragueSkyModel.cpp
    PragueSkyModelTest.cpp
)

add_executable(PragueSkyModel ${PRAGUESKYMODEL_SRC})