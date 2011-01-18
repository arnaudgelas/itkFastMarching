project( FMM )

cmake_minimum_required( 2.8 )

find_package( ITK REQUIRED )

if( ${ITK_FOUND} )
  include( ${ITK_USE_FILE} )
endif()

include_directories( Code )

add_subdirectory( Testing )

