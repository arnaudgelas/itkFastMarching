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

# --------------------------------------------------------------------------Topology constrained front propagation
# --------------------------------------------------------------------------

# non-templated class -- this should be stored in a library and linked in...
ADD_EXECUTABLE( itkFastMarchingImageTopologicalTest
  itkFastMarchingImageTopologicalTest.cxx
  ${BASE_SOURCE}
)
TARGET_LINK_LIBRARIES( itkFastMarchingImageTopologicalTest ITKCommon ITKIO )

## 2-D

ADD_TEST( TEST_Brain2D_singleSeed_NoTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 2
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_speed.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_Brain2D_singleSeed_NoTopo.nii.gz
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_singleSeed.nii.gz
  150 0 )

ADD_TEST( TEST_Brain2D_singleSeed_StrictTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 2
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_speed.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_Brain2D_singleSeed_StrictTopo.nii.gz
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_singleSeed.nii.gz
  150 1 )

ADD_TEST( TEST_Brain2D_singleSeed_NoHandlesTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 2
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_speed.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_Brain2D_singleSeed_NoHandlesTopo.nii.gz
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_singleSeed.nii.gz
  150 2 )

ADD_TEST( TEST_Brain2D_multipleSeeds_NoTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 2
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_speed.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_Brain2D_multipleSeeds_NoTopo.nii.gz
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_multipleSeeds.nii.gz
  150 0 )

ADD_TEST( TEST_Brain2D_multipleSeeds_StrictTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 2
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_speed.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_Brain2D_multipleSeeds_StrictTopo.nii.gz
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_multipleSeeds.nii.gz
  150 1 )

ADD_TEST( TEST_Brain2D_multipleSeeds_NoHandlesTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 2
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_speed.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_Brain2D_multipleSeeds_NoHandlesTopo.nii.gz
  ${DATA_DIRECTORY}/BrainProtonDensitySlice_multipleSeeds.nii.gz
  150 2 )

## 3-D

ADD_TEST( TEST_torus_multipleSeeds_NoTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 3
  ${DATA_DIRECTORY}/torus.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_torus_multipleSeeds_NoTopo.nii.gz
  ${DATA_DIRECTORY}/torus_multipleSeeds.nii.gz
  150 0 )

ADD_TEST( TEST_torus_multipleSeeds_StrictTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 3
  ${DATA_DIRECTORY}/torus.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_torus_multipleSeeds_StrictTopo.nii.gz
  ${DATA_DIRECTORY}/torus_multipleSeeds.nii.gz
  150 1 )

ADD_TEST( TEST_torus_multipleSeeds_NoHandlesTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 3
  ${DATA_DIRECTORY}/torus.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_torus_multipleSeeds_NoHandlesTopo.nii.gz
  ${DATA_DIRECTORY}/torus_multipleSeeds.nii.gz
  150 2 )

ADD_TEST( TEST_wm_multipleSeeds_NoTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 3
  ${DATA_DIRECTORY}/wm.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_wm_multipleSeeds_NoTopo.nii.gz
  ${DATA_DIRECTORY}/wm_multipleSeeds.nii.gz
  150 0 )

ADD_TEST( TEST_wm_multipleSeeds_StrictTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 3
  ${DATA_DIRECTORY}/wm.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_wm_multipleSeeds_StrictTopo.nii.gz
  ${DATA_DIRECTORY}/wm_multipleSeeds.nii.gz
  150 1 )

ADD_TEST( TEST_wm_multipleSeeds_NoHandlesTopo
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFastMarchingImageFilterTest 3
  ${DATA_DIRECTORY}/wm.nii.gz
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_wm_multipleSeeds_NoHandlesTopo.nii.gz
  ${DATA_DIRECTORY}/wm_multipleSeeds.nii.gz
  150 2 )

