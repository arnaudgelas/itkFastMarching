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

#include "itkFastMarchingImageFilterBase.h"
#include "itkFastMarchingThresholdStoppingCriterion.h"
#include "itkImage.h"

int main( int argc, char** argv )
  {
  (void) argc;
  (void) argv;

  const unsigned int Dimension2 = 2;
  typedef float PixelType;

  typedef itk::Image< PixelType, Dimension2 > Image2DType;
  Image2DType::Pointer input2d = Image2DType::New();

  typedef itk::FastMarchingImageThresholdStoppingCriterion< Image2DType >
      CriterionType2D;

  typedef itk::FastMarchingImageFilterBase<
      Dimension2,
      PixelType,
      PixelType,
      CriterionType2D > FMM2DType;
  FMM2DType::Pointer fmm2d = FMM2DType::New();
  fmm2d->SetInput( input2d );
  fmm2d->Update();

  Image2DType::Pointer output2d = fmm2d->GetOutput();

  const unsigned int Dimension3 = 3;
  typedef itk::Image< float, Dimension3 > Image3DType;
  Image3DType::Pointer input3d = Image3DType::New();

  typedef itk::FastMarchingImageThresholdStoppingCriterion< Image3DType >
      CriterionType3D;

  typedef itk::FastMarchingImageFilterBase< Dimension3,
      float,
      float,
      CriterionType3D > FMM3DType;
  FMM3DType::Pointer fmm3d = FMM3DType::New();
  fmm3d->SetInput( input3d );
  fmm3d->Update();

  Image3DType::Pointer output3d = fmm3d->GetOutput();

  return EXIT_SUCCESS;
  }
