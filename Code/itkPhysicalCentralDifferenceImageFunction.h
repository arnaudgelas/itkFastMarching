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

#ifndef __itkPhysicalCentralDifferenceImageFunction_h
#define __itkPhysicalCentralDifferenceImageFunction_h

#include "itkImageFunction.h"
#include "itkCovariantVector.h"
#include "itkImageBase.h"
#include "itkLinearInterpolateImageFunction.h"

namespace itk
{

/**
 * \class PhysicalCentralDifferenceImageFunction
 * \brief Calculate the derivative by central differencing in physical space.
 *
 * This class is templated over the input image type and
 * the coordinate representation type (e.g. float or double).
 *
 * \author Dan Mueller, Queensland University of Technology, dan.muel[at]gmail.com
 *
 * \ingroup ImageFunctions
 */
template <
  class TInputImage,
  class TCoordRep = float >
class ITK_EXPORT PhysicalCentralDifferenceImageFunction :
  public ImageFunction< TInputImage,
                        CovariantVector<double, \
                        ::itk::GetImageDimension<TInputImage>::ImageDimension>,
                        TCoordRep >
{
public:
  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard class typedefs. */
  typedef PhysicalCentralDifferenceImageFunction Self;
  typedef ImageFunction<TInputImage,
                        CovariantVector<double,
                        itkGetStaticConstMacro(ImageDimension)>,
                        TCoordRep>       Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(PhysicalCentralDifferenceImageFunction, ImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** InputImageType typedef support. */
  typedef TInputImage InputImageType;

  /** OutputType typdef support. */
  typedef typename Superclass::OutputType OutputType;

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Point typedef support. */
  typedef typename Superclass::PointType PointType;

  /** Linear interpolate function typedef. */
  typedef LinearInterpolateImageFunction<TInputImage, TCoordRep>
     InterpolateImageFunctionType;

  /** Set the input image.
   * \warning this method caches BufferedRegion information.
   * If the BufferedRegion has changed, user must call
   * SetInputImage again to update cached values. */
  virtual void SetInputImage( const InputImageType * ptr )
    {
    this->Superclass::SetInputImage( ptr );
    if ( m_Interpolator.IsNotNull() )
      {
      m_Interpolator->SetInputImage( ptr );
      }
    }

  /** Evalulate the image derivative by central differencing at specified index.
   *
   *  No bounds checking is done. The point is assume to lie within the
   *  image buffer. ImageFunction::IsInsideBuffer() can be used to check
   *  bounds before calling this method. */
  virtual OutputType EvaluateAtIndex( const IndexType& index ) const
    {
    PointType point;
    m_Interpolator->GetInputImage()->TransformIndexToPhysicalPoint( index, point );
    return this->Evaluate( point );
    }

  /** Evalulate the image derivative by central differencing at specified
   *  continuous index.
   *
   *  No bounds checking is done. The point is assume to lie within the
   *  image buffer. ImageFunction::IsInsideBuffer() can be used to check
   *  bounds before calling this method. */
  virtual OutputType EvaluateAtContinuousIndex(
    const ContinuousIndexType& cindex ) const
    {
    PointType point;
    m_Interpolator->GetInputImage()->TransformContinuousIndexToPhysicalPoint(
          cindex, point );
    return this->Evaluate( point );
    }

  /** Evalulate the image derivative by central differencing at specified
   *  physical point.
   *
   *  No bounds checking is done. The point is assume to lie within the
   *  image buffer. ImageFunction::IsInsideBuffer() can be used to check
   *  bounds before calling this method. */
   virtual OutputType Evaluate( const PointType& point ) const;

protected:
  PhysicalCentralDifferenceImageFunction();
  ~PhysicalCentralDifferenceImageFunction(){};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  PhysicalCentralDifferenceImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  typename InterpolateImageFunctionType::Pointer m_Interpolator;

};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_PhysicalCentralDifferenceImageFunction(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT PhysicalCentralDifferenceImageFunction< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef PhysicalCentralDifferenceImageFunction< ITK_TEMPLATE_2 x > \
                        PhysicalCentralDifferenceImageFunction##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkPhysicalCentralDifferenceImageFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkPhysicalCentralDifferenceImageFunction.txx"
#endif

#endif
