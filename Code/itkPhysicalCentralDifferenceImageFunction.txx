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

#ifndef __itkPhysicalCentralDifferenceImageFunction_txx
#define __itkPhysicalCentralDifferenceImageFunction_txx

#include "itkPhysicalCentralDifferenceImageFunction.h"

namespace itk
{

template <class TInputImage, class TCoordRep>
PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>
::PhysicalCentralDifferenceImageFunction()
{
  m_Interpolator = InterpolateImageFunctionType::New();
}

template <class TInputImage, class TCoordRep>
void
PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
}

template <class TInputImage, class TCoordRep>
typename PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>::OutputType
PhysicalCentralDifferenceImageFunction<TInputImage,TCoordRep>
::Evaluate( const PointType& point ) const
{
  OutputType derivative;
  derivative.Fill( 0.0 );

  typename InputImageType::SpacingType spacing = this->m_Image->GetSpacing();
  double s;

  for(unsigned int dim=0; dim < ImageDimension; dim++)
    {
    s = static_cast< double >( spacing[dim] );

    // Get the left neighbor
    PointType pointLeft( point );
    pointLeft[dim] -= s;
    TCoordRep valueLeft = m_Interpolator->Evaluate( pointLeft );

    // Get the right neighbor
    PointType pointRight( point );
    pointRight[dim] += s;
    TCoordRep valueRight = m_Interpolator->Evaluate( pointRight );

    // Compute derivative
    derivative[dim] = (valueRight - valueLeft) / ( 2. * s );
    }

  return ( derivative );
}

} // end namespace itk

#endif
