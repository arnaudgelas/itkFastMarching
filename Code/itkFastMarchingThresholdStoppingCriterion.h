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

#ifndef __itkFastMarchingThresholdStoppingCriterion_h
#define __itkFastMarchingThresholdStoppingCriterion_h

#include "itkFastMarchingStoppingCriterionBase.h"
#include "itkObjectFactory.h"

namespace itk
{
  /**
    \class FastMarchingThresholdStoppingCriterion
    \brief Stopping Criterion is verified when Current Value is Above
    the provided threshold.
    */
  template< class TNode,
           typename TValue >
  class FastMarchingThresholdStoppingCriterion :
      public FastMarchingStoppingCriterionBase< TNode, TValue>
    {
  public:
    typedef FastMarchingThresholdStoppingCriterion Self;
    typedef FastMarchingStoppingCriterionBase< TNode, TValue > Superclass;
    typedef SmartPointer< Self >              Pointer;
    typedef SmartPointer< const Self >        ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(FastMarchingThresholdStoppingCriterion,
                 FastMarchingStoppingCriterionBase );

    typedef typename Superclass::ValueType  ValueType;
    typedef typename Superclass::NodeType   NodeType;

    itkSetMacro( Threshold, ValueType );
    itkGetMacro( Threshold, ValueType );

    void SetCurrentNode( const NodeType& iNode ) { (void) iNode; }
    bool IsSatisfied() const
      {
      return ( this->m_CurrentValue >= this->m_Threshold );
      }

    const std::string GetDescription() const
      {
      return "Current Value >= Threshold";
      }

  protected:
    FastMarchingThresholdStoppingCriterion() : Superclass()
    {
      m_Threshold = NumericTraits< ValueType >::Zero;
    }
    ~FastMarchingThresholdStoppingCriterion() {}

    ValueType m_Threshold;

  private:
    FastMarchingThresholdStoppingCriterion( const Self& );
    void operator = ( const Self& );
    };

  /** \class FastMarchingImageThresholdStoppingCriterion
    \brief Partial Specialization of FastMarchingThresholdStoppingCriterion
    for Image */
  template< class TImage >
  class FastMarchingImageThresholdStoppingCriterion :
    public FastMarchingThresholdStoppingCriterion<
      typename TImage::IndexType,
      typename TImage::PixelType >
    {
  public:
    typedef TImage ImageType;
    typedef typename ImageType::IndexType NodeType;
    typedef typename ImageType::PixelType ValueType;

    typedef FastMarchingImageThresholdStoppingCriterion Self;
    typedef FastMarchingThresholdStoppingCriterion< NodeType,
      ValueType > Superclass;
    typedef SmartPointer< Self >              Pointer;
    typedef SmartPointer< const Self >        ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(FastMarchingImageThresholdStoppingCriterion,
                 FastMarchingThresholdStoppingCriterion );

  protected:
    FastMarchingImageThresholdStoppingCriterion() : Superclass() {}
    ~FastMarchingImageThresholdStoppingCriterion() {}
    };

}
#endif // __itkFastMarchingThresholdStoppingCriterion_h
