project( FMM )

cmake_minimum_required(VERSION 2.8)

find_package( ITK REQUIRED )

if( ${ITK_FOUND} )
  include( ${ITK_USE_FILE} )
endif()

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${FMM_BINARY_DIR}/bin
  )

set( DATA_DIRECTORY
  ${FMM_SOURCE_DIR}/Data
  )

set( CTEST_CUSTOM_COVERAGE_EXCLUDE
  ${CTEST_CUSTOM_COVERAGE_EXCLUDE}
  "/Testing/*"
)

set( BASE_SOURCE ${FMM_SOURCE_DIR}/Code/itkStoppingCriterionBase.cxx )

include_directories( Code )

include( CTest )
add_subdirectory( Testing )

if( BUILD_TESTING )
  set(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of the build on the dashboard" )
endif()
