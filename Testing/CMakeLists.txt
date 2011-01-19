add_executable( itkFastMarchingBaseTest itkFastMarchingBaseTest.cxx )
target_link_libraries( itkFastMarchingBaseTest ITKCommon )

add_test( itkImageFastMarchingBaseTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingBaseTest 0 )

add_test( itkMeshFastMarchingBaseTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingBaseTest 1 )

add_executable( itkFastMarchingImageFilterBaseTest
 itkFastMarchingImageFilterBaseTest.cxx )
target_link_libraries( itkFastMarchingImageFilterBaseTest ITKCommon )
