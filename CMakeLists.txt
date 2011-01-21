project( FMM )

cmake_minimum_required(VERSION 2.8)

find_package( ITK REQUIRED )

if( ${ITK_FOUND} )
  include( ${ITK_USE_FILE} )
endif()

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${FMM_BINARY_DIR}/bin
  )

set( BASE_SOURCE ${FMM_SOURCE_DIR}/Code/itkStoppingCriterionBase.cxx )

message( ${BASE_SOURCE} )

include_directories( Code )

include( CTest )
add_subdirectory( Testing )

