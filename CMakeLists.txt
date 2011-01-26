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
  "${FMM_SOURCE_DIR}/Testing/*"
)

set( BASE_SOURCE ${FMM_SOURCE_DIR}/Code/itkStoppingCriterionBase.cxx )

include_directories( Code )

option( BUILD_DOCUMENTATION "build doxygen doc" OFF )

if( BUILD_DOCUMENTATION )
  find_package( Doxygen )

  configure_file( ${FMM_SOURCE_DIR}/doxygen.config.in
    ${FMM_BINARY_DIR}/doxygen.config)


  add_custom_target(Documentation ALL
    COMMAND ${DOXYGEN_EXECUTABLE} ${FMM_BINARY_DIR}/doxygen.config
    MAIN_DEPENDENCY ${FMM_BINARY_DIR}/doxygen.config
    WORKING_DIRECTORY ${FMM_BINARY_DIR}/
  )
endif()

include( CTest )
add_subdirectory( Testing )

if( BUILD_TESTING )
  set(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of the build on the dashboard" )
endif()
