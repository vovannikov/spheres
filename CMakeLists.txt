CMAKE_MINIMUM_REQUIRED(VERSION 3.13.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(TARGET "bone")

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

FIND_PACKAGE(deal.II 9.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )

IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** BONE: Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(${TARGET} VERSION 1.0 LANGUAGES CXX)

include_directories(${DEAL_II_DIR}/include)

#include_directories(/usr/lib/x86_64-linux-gnu/openmpi/include)

include(FetchContent)
if(${CMAKE_VERSION} VERSION_LESS 3.14)
    include(cmake/add_FetchContent_MakeAvailable.cmake)
endif()

FetchContent_Declare(
  pfsim
  GIT_REPOSITORY git@github.com:vovannikov/pfsim.git
  GIT_TAG        pf_monolithic
  SOURCE_SUBDIR  sources
)
FetchContent_MakeAvailable(pfsim)

add_subdirectory(src)
