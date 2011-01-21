/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkFastMarchingThresholdStoppingCriterion.h"
#include "itkFastMarchingBase.h"

int main( int argc, char* argv[] )
  {
  (void) argc;
  (void) argv;

   typedef itk::ImageFastMarchingTraits< 2, float, float >
    ImageTraits;
  typedef itk::FastMarchingThresholdStoppingCriterion< ImageTraits >
    ImageStoppingCriterionType;

  ImageStoppingCriterionType::Pointer image_criterion =
      ImageStoppingCriterionType::New();

  typedef itk::MeshFastMarchingTraits< 3,
      float,
      itk::DefaultStaticMeshTraits< float, 3, 3 >,
      float,
      itk::DefaultStaticMeshTraits< float, 3, 3 > >
    MeshTraits;

  typedef itk::FastMarchingThresholdStoppingCriterion< MeshTraits >
    MeshStoppingCriterionType;

  MeshStoppingCriterionType::Pointer mesh_criterion =
      MeshStoppingCriterionType::New();

  return EXIT_SUCCESS;
  }
