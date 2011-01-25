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
#ifndef __itkFastMarchingExtensionImageFilterBase_h
#define __itkFastMarchingExtensionImageFilterBase_h

#include "itkFastMarchingImageFilterBase.h"

namespace itk
{
/** \class FastMarchingExtensionImageFilterBase
 * \brief Extend auxiliary variables smoothly using Fast Marching.
 *
 * Fast marching can be used to extend auxiliary variables smoothly
 * from the zero level set. Starting from an initial position on the
 * front, this class simultaneously calculate the signed distance and
 * extend a set of auxiliary values.
 *
 * This class is templated over the level set image type, the auxiliary
 * variable type and the number of auxiliary variables to extended. The initial
 * front is specified by two containers: one containing the known points
 * and one containing the trial points. The auxiliary variables on the front
 * are represented by two auxiliary variable containers: one containing
 * the value of the variables at the know points and on containing the
 * value of the variables at the trail points.
 *
 * Implemenation of this class is based on Chapter 11 of
 * "Level Set Methods and Fast Marching Methods", J.A. Sethian,
 * Cambridge Press, Second edition, 1999.
 *
 * \sa FastMarchingImageFilter
 * \sa LevelSetTypeDefault
 * \sa AuxVarTypeDefault
 * \ingroup LevelSetSegmentation
 */
template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TCriterion,
         typename TAuxValue,
         unsigned int VAuxDimension >
class ITK_EXPORT FastMarchingExtensionImageFilterBase:
  public FastMarchingImageFilterBase< VDimension, TInputPixel,
    TOutputPixel, TCriterion >
{
public:
  /** Standard class typdedefs. */
  typedef FastMarchingExtensionImageFilterBase                  Self;
  typedef FastMarchingImageFilterBase< VDimension, TInputPixel,
    TOutputPixel, TCriterion > Superclass;
  typedef SmartPointer< Self >                              Pointer;
  typedef SmartPointer< const Self >                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FastMarchingExtensionImageFilterBase,
               FastMarchingImageFilterBase);

  /** The dimension of the level set. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension );

  /** Number of auxiliary variables to be extended. */
  itkStaticConstMacro(AuxDimension, unsigned int, VAuxDimension);

  /** AuxVarType typedef support. */
  //typedef AuxVarTypeDefault< TAuxValue,
  //                           itkGetStaticConstMacro(AuxDimension),
  //                           itkGetStaticConstMacro(ImageDimension) >
  //AuxVarType;
  typedef TAuxValue       AuxValueType;
  typedef Vector< AuxValueType, AuxDimension >  AuxValueVectorType;
  typedef std::vector< AuxValueVectorType >     AuxValueContainer;
  typedef typename AuxValueContainer::const_iterator AuxValueContainerConstIterator;

  typedef Image< AuxValueType, ImageDimension > AuxImageType;
  typedef typename AuxImageType::Pointer        AuxImagePointer;

  //typedef typename AuxVarType::AuxValueContainer  AuxValueContainer;
  //typedef typename AuxVarType::AuxImageType       AuxImageType;
  //typedef typename AuxVarType::AuxImagePointer    AuxImagePointer;

  /** Index typedef support. */
  typedef typename Superclass::NodeType NodeType;
  typedef typename Superclass::NodePairType NodePairType;
  typedef typename Superclass::NodeContainerType NodeContainerType;
  typedef typename Superclass::NodeContainerConstIterator NodeContainerConstIterator;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename Superclass::OutputPixelType OutputPixelType;
  typedef typename Superclass::InternalNodeStructure InternalNodeStructure;

  /** Get one of the extended auxiliary variable image. */
  AuxImageType * GetAuxiliaryImage(unsigned int idx);

  /** Set the container auxiliary values at the initial alive points. */
  void SetAuxiliaryAliveValues(AuxValueContainer values)
  {
    m_AuxAliveValues = values;
  }

  /** Get the container of auxiliary values at the initial alive points. */
  AuxValueContainer GetAuxiliaryAliveValues(void)
  {
    return m_AuxAliveValues;
  }

  /** Set the container of auxiliary values at the initial trial points. */
  void SetAuxiliaryTrialValues(AuxValueContainer values)
  {
    m_AuxTrialValues = values;
  }

  /** Get the container of auxiliary values at the initial trial points. */
  AuxValueContainer GetAuxiliaryTrialValues()
  {
    return m_AuxTrialValues;
  }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( AuxValueHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< TAuxValue > ) );
  /** End concept checking */
#endif
protected:
  FastMarchingExtensionImageFilterBase();
  ~FastMarchingExtensionImageFilterBase(){}
  void PrintSelf(std::ostream & os, Indent indent) const;

  virtual void InitializeOutput(OutputImageType *);

  virtual void UpdateValue( OutputImageType* oImage, const NodeType& iValue );

  /** Generate the output image meta information */
  virtual void GenerateOutputInformation();

  virtual void EnlargeOutputRequestedRegion(DataObject *output);

private:
  FastMarchingExtensionImageFilterBase(const Self &); //purposely not implemented
  void operator=(const Self &);                   //purposely not implemented

  AuxValueContainer m_AuxAliveValues;
  AuxValueContainer m_AuxTrialValues;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFastMarchingExtensionImageFilterBase.txx"
#endif

#endif
