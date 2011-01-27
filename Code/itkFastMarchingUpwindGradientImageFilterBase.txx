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
#ifndef __itkFastMarchingUpwindGradientImageFilterBase_txx
#define __itkFastMarchingUpwindGradientImageFilterBase_txx

#include "itkFastMarchingUpwindGradientImageFilterBase.h"
#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include <algorithm>

namespace itk
{
/**
 *
 */
template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TSuperclass >
FastMarchingUpwindGradientImageFilterBase< VDimension, TInputPixel, TOutputPixel, TSuperclass >
::FastMarchingUpwindGradientImageFilterBase()
{
  m_GradientImage = GradientImageType::New();
  m_GenerateGradientImage = false;
}

/**
 *
 */
template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TSuperclass >
void
FastMarchingUpwindGradientImageFilterBase< VDimension, TInputPixel, TOutputPixel, TSuperclass >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Gradient image: " << m_GradientImage.GetPointer() << std::endl;
  os << indent << "Generate gradient image: " << m_GenerateGradientImage << std::endl;
}

/**
 *
 */
template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TSuperclass >
void
FastMarchingUpwindGradientImageFilterBase< VDimension, TInputPixel, TOutputPixel, TSuperclass >::
InitializeOutput(OutputImageType *output)
{
  Superclass::InitializeOutput(output);

  // allocate memory for the GradientImage if requested
  if ( m_GenerateGradientImage )
    {
    m_GradientImage->CopyInformation( this->GetInput() );
    m_GradientImage->SetBufferedRegion( output->GetBufferedRegion() );
    m_GradientImage->Allocate();
    }

  // set all gradient vectors to zero
  if ( m_GenerateGradientImage )
    {
    typedef ImageRegionIterator< GradientImageType > GradientIterator;

    GradientIterator gradientIt( m_GradientImage,
                                 m_GradientImage->GetBufferedRegion() );

    GradientPixelType zeroGradient;
    typedef typename GradientPixelType::ValueType GradientPixelValueType;
    zeroGradient.Fill(NumericTraits< GradientPixelValueType >::Zero);

    gradientIt.GoToBegin();

    while( !gradientIt.IsAtEnd() )
      {
      gradientIt.Set(zeroGradient);
      ++gradientIt;
      }
    }
}

template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TSuperclass >
void
FastMarchingUpwindGradientImageFilterBase< VDimension, TInputPixel, TOutputPixel, TSuperclass >::
UpdateNeighbors(
  OutputImageType* oImage,
  const NodeType& iNode )
{
  Superclass::UpdateNeighbors( oImage, iNode );

  if ( m_GenerateGradientImage )
    {
    this->ComputeGradient( oImage, iNode );
    }
}

/**
 *
 */
template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TSuperclass >
void
FastMarchingUpwindGradientImageFilterBase< VDimension, TInputPixel, TOutputPixel, TSuperclass >
::ComputeGradient( OutputImageType* oImage,
                  const NodeType& iNode )
{
  NodeType neighIndex = iNode;

  OutputPixelType centerPixel;
  OutputPixelType dx_forward;
  OutputPixelType dx_backward;
  GradientPixelType gradientPixel;

  const OutputPixelType ZERO = NumericTraits< OutputPixelType >::Zero;

  OutputSpacingType spacing = oImage->GetSpacing();

  unsigned int xStride[ itkGetStaticConstMacro( ImageDimension ) ];

  centerPixel = oImage->GetPixel( iNode );

  for ( unsigned int j = 0; j < ImageDimension; j++ )
    {
    neighIndex = iNode;

    // Set stride of one in each direction
    xStride[j] = 1;

    // Compute one-sided finite differences with alive neighbors
    // (the front can only come from there)
    dx_backward = ZERO;
    neighIndex[j] = iNode[j] - xStride[j];

    if ( !( ( neighIndex[j] > this->m_LastIndex[j] ) ||
            ( neighIndex[j] < this->m_StartIndex[j] ) ) )
      {
      if ( this->GetLabelValueForGivenNode(neighIndex) == Superclass::Alive )
        {
        dx_backward = centerPixel - oImage->GetPixel(neighIndex);
        }
      }

    dx_forward = ZERO;
    neighIndex[j] = iNode[j] + xStride[j];

    if ( !( ( neighIndex[j] > this->m_LastIndex[j] ) ||
            ( neighIndex[j] < this->m_StartIndex[j] ) ) )
      {
      if ( this->GetLabelValueForGivenNode(neighIndex) == Superclass::Alive )
        {
        dx_forward = oImage->GetPixel(neighIndex) - centerPixel;
        }
      }

    // Compute upwind finite differences
    if ( vnl_math_max(dx_backward, -dx_forward) < ZERO )
      {
      gradientPixel[j] = ZERO;
      }
    else
      {
      if ( dx_backward > -dx_forward )
        {
        gradientPixel[j] = dx_backward;
        }
      else
        {
        gradientPixel[j] = dx_forward;
        }
      }

    gradientPixel[j] /= spacing[j];
    }

  m_GradientImage->SetPixel(iNode, gradientPixel);
}
} // namespace itk

#endif
