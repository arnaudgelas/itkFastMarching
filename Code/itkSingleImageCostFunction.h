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

#ifndef __itkSingleImageCostFunction_h
#define __itkSingleImageCostFunction_h

#include "itkNumericTraits.h"
#include "itkExceptionObject.h"
#include "itkContinuousIndex.h"
#include "itkSingleValuedCostFunction.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkCentralDifferenceImageFunction.h"
#include "itkPhysicalCentralDifferenceImageFunction.h"

namespace itk
{

/** \class SingleImageCostFunction
 * \brief This class is a cost function which queries
 *        an underlying image for the single value.
 *
 * The user is expected to provide an image representing the
 * underlying cost function. The user may also provide an image
 * interpolator (if not provided the LinearInterpolateImageFunction
 * is used by default). The gradient is computed using central
 * differences in physical space.
 *
 * The parameters are the physical location (itkPoint) of the current
 * position. Initialize() must be called before using this cost function.
 *
 * \author Dan Mueller, Queensland University of Technology, dan.muel[at]gmail.com
 *
 * \ingroup Numerics Optimizers
 */
template <class TImage>
class ITK_EXPORT SingleImageCostFunction :
    public SingleValuedCostFunction
{
public:
  /** Standard class typedefs. */
  typedef SingleImageCostFunction      Self;
  typedef SingleValuedCostFunction     Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( SingleImageCostFunction, SingleValuedCostFunction );

  /** MeasureType typedef.
   *  It defines a type used to return the cost function value. */
  typedef typename Superclass::MeasureType          MeasureType;

  /** DerivativeType typedef.
   *  It defines a type used to return the cost function derivative. */
  typedef typename Superclass::DerivativeType       DerivativeType;

  /** ParametersType typedef.
   *  It defines a position in the optimization search space. */
  typedef typename Superclass::ParametersType       ParametersType;

  /**  Type of the Image. */
  typedef TImage                           ImageType;
  typedef typename TImage::PixelType       ImagePixelType;
  typedef typename ImageType::ConstPointer ImageConstPointer;

  /** Constant for the image dimension */
  itkStaticConstMacro(ImageDimension, unsigned int, ImageType::ImageDimension);

  /** Type used for representing point components */
  typedef Superclass::ParametersValueType CoordRepType;

  /** Type for locations */
  typedef Index< ImageDimension > IndexType;
  typedef Point< CoordRepType, ImageDimension > PointType;
  typedef ContinuousIndex< CoordRepType, ImageDimension >
    ContinuousIndexType;

  /** Type of the Interpolator class */
  typedef InterpolateImageFunction< ImageType, CoordRepType >
    InterpolatorType;
  typedef LinearInterpolateImageFunction< ImageType, CoordRepType >
    DefaultInterpolatorType;

  /** Type of the GradientImageFunction class */
  typedef PhysicalCentralDifferenceImageFunction< ImageType, CoordRepType >
    GradientImageFunctionType;

  /** Get/set the Interpolator. */
  itkSetObjectMacro( Interpolator, InterpolatorType );
  itkGetConstObjectMacro( Interpolator, InterpolatorType );

  /** Get/set the Image.  */
  itkSetConstObjectMacro( Image, ImageType );
  itkGetConstObjectMacro( Image, ImageType );

  /** Initialize the cost function */
  virtual void Initialize(void) throw ( ExceptionObject );

  /** Return the number of parameters required by the Transform */
  unsigned int GetNumberOfParameters(void) const
  { return ImageDimension; }

  /** This method returns the value of the cost function corresponding
    * to the specified parameters. */
  virtual MeasureType GetValue( const ParametersType & parameters ) const;

  /** This method returns the derivative of the cost function corresponding
    * to the specified parameters. */
  virtual void GetDerivative( const ParametersType & parameters,
                              DerivativeType & derivative ) const;

protected:
  SingleImageCostFunction();
  virtual ~SingleImageCostFunction() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  SingleImageCostFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  ImageConstPointer                           m_Image;
  typename InterpolatorType::Pointer          m_Interpolator;
  typename GradientImageFunctionType::Pointer m_GradientImageFunction;
  typename DerivativeType::ValueType          m_DerivativeThreshold;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSingleImageCostFunction.txx"
#endif

#endif
