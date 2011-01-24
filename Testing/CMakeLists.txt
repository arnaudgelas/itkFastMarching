add_executable( itkFastMarchingBaseTest
  itkFastMarchingBaseTest.cxx
)
target_link_libraries( itkFastMarchingBaseTest ITKCommon )

add_test( itkImageFastMarchingBaseTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingBaseTest 0 )

add_test( itkMeshFastMarchingBaseTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingBaseTest 1 )

add_executable( itkFastMarchingImageFilterBaseTest
  itkFastMarchingImageFilterBaseTest.cxx
)
target_link_libraries( itkFastMarchingImageFilterBaseTest ITKCommon )

add_test( itkFastMarchingImageFilterBaseTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterBaseTest )


add_executable( itkFastMarchingStoppingCriterionBaseTest
  itkFastMarchingStoppingCriterionBaseTest.cxx
  ${BASE_SOURCE}
)
target_link_libraries( itkFastMarchingStoppingCriterionBaseTest ITKCommon )


add_executable( itkFastMarchingThresholdStoppingCriterionTest
  itkFastMarchingThresholdStoppingCriterionTest.cxx
  ${BASE_SOURCE}
)
target_link_libraries( itkFastMarchingThresholdStoppingCriterionTest ITKCommon )

add_executable( itkFastMarchingImageFilterRealTest1
  itkFastMarchingImageFilterRealTest1.cxx
  ${BASE_SOURCE}
)
target_link_libraries( itkFastMarchingImageFilterRealTest1 ITKCommon )


add_test( itkFastMarchingImageFilterRealTest1
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterRealTest1 )
