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
  template< class TTraits >
  class FastMarchingThresholdStoppingCriterion :
      public FastMarchingStoppingCriterionBase< TTraits >
    {
  public:
    typedef FastMarchingThresholdStoppingCriterion Self;
    typedef FastMarchingStoppingCriterionBase< TTraits > Superclass;
    typedef SmartPointer< Self >              Pointer;
    typedef SmartPointer< const Self >        ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(FastMarchingThresholdStoppingCriterion,
                 FastMarchingStoppingCriterionBase );

    typedef typename Superclass::OutputPixelType  OutputPixelType;
    typedef typename Superclass::NodeType         NodeType;

    itkSetMacro( Threshold, OutputPixelType );
    itkGetMacro( Threshold, OutputPixelType );

    void SetCurrentNode( const NodeType& iNode ) { (void) iNode; }
    bool IsSatisfied() const
      {
      return ( this->m_CurrentValue >= this->m_Threshold );
      }

  protected:
    FastMarchingThresholdStoppingCriterion() : Superclass()
    {
      m_Threshold = NumericTraits< OutputPixelType >::Zero;
    }
    ~FastMarchingThresholdStoppingCriterion() {}

    OutputPixelType m_Threshold;

  private:
    FastMarchingThresholdStoppingCriterion( const Self& );
    void operator = ( const Self& );
    };
}
#endif // __itkFastMarchingThresholdStoppingCriterion_h
