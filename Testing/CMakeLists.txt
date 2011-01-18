add_executable( itkFastMarchingBaseTest itkFastMarchingBaseTest.cxx )
target_link_libraries( itkFastMarchingBaseTest ITKCommon )

add_test( itkFastMarchingBaseTest itkFastMarchingBaseTest )
