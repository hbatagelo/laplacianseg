#############################################################################
##
## Copyright (C) 2020-2022 Harlen Batagelo <hbatagelo@gmail.com>,
##                         João Paulo Gois <jpgois@gmail.com>.
##
## This file is part of the implementation of the paper
## 'Laplacian Coordinates: Theory and Methods for Seeded Image Segmentation'
## by Wallace Casaca, João Paulo Gois, Harlen Batagelo, Gabriel Taubin and
## Luis Gustavo Nonato. DOI 10.1109/TPAMI.2020.2974475.
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.
##
#############################################################################

cmake_minimum_required(VERSION 3.16)

project(lcseg VERSION 1.0.1 LANGUAGES CXX)

set(CMAKE_AUTOMOC ON)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # Suppress Eigen's anonymous enum warnings See
  # https://gitlab.com/libeigen/eigen/-/issues/1826
  add_compile_options(-Wno-deprecated-enum-enum-conversion)
endif()

find_package(Qt6 COMPONENTS Core Concurrent Gui REQUIRED)
add_library(Qt STATIC)
target_link_libraries(Qt Qt6::Core Qt6::Concurrent Qt6::Gui)

set(SOURCES commandlineparser.cpp main.cpp seededsegmentationhard.cpp
            seededsegmentationsoft.cpp seededspsegmentation.cpp slic.cpp)
set(HEADERS commandlineparser.h seededsegmentationhard.h
            seededsegmentationsoft.h seededspsegmentation.h slic.h util.h)

find_package(Eigen3 3.3 REQUIRED)
find_package(OpenCV 4.2 COMPONENTS core imgproc REQUIRED)
find_package(SuiteSparse REQUIRED)
find_package(OpenMP)

if(SuiteSparse_CHOLMOD_FOUND)
  message(STATUS "Found SuiteSparse CHOLMOD: "
          ${SuiteSparse_CHOLMOD_LIBRARY_RELEASE})
endif()

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

target_include_directories(${PROJECT_NAME} SYSTEM
                           PRIVATE ${EIGEN3_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} SYSTEM
                           PRIVATE ${OpenCV_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} SYSTEM
                           PRIVATE ${SuiteSparse_INCLUDE_DIRS})

target_link_libraries(
  ${PROJECT_NAME} PRIVATE Qt ${OpenCV_LIBS} ${SuiteSparse_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)
