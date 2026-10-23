cmake_minimum_required(VERSION 3.16)

# Portable DSP sources are shared with the host test build.
set(EXTRA_COMPONENT_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../dsp")

include($ENV{IDF_PATH}/tools/cmake/project.cmake)
project(ecg_monitor)
