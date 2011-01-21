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

#ifndef __itkFastMarchingStoppingCriterionBase_h
#define __itkFastMarchingStoppingCriterionBase_h

#include "itkStoppingCriterionBase.h"
#include "itkNumericTraits.h"

namespace itk
{

 /** \class FastMarchingStoppingCriterionBase

 */
 template< class TTraits >
 class FastMarchingStoppingCriterionBase : public StoppingCriterionBase
 {
 public:
  typedef FastMarchingStoppingCriterionBase Self;
  typedef StoppingCriterionBase             Superclass;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;

  typedef TTraits Traits;
  typedef typename Traits::OutputPixelType  OutputPixelType;
  typedef typename Traits::NodeType         NodeType;



  /** Run-time type information (and related methods). */
  itkTypeMacro(FastMarchingStoppingCriterionBase, StoppingCriterionBase);

  virtual void SetCurrentNode( const NodeType& iNode ) = 0;

  void SetCurrentValue( const OutputPixelType& iValue )
  {
    if( iValue >= m_CurrentValue )
      {
      m_PreviousValue = m_CurrentValue;
      m_CurrentValue = iValue;
      }
    else
      {
      itkExceptionMacro( << "Current value is decreasing!" );
      }
  }

 protected:
  FastMarchingStoppingCriterionBase() : Superclass()
  {

    m_CurrentValue = NumericTraits< OutputPixelType >::Zero;
    m_PreviousValue = NumericTraits< OutputPixelType >::Zero;
  }
  virtual ~FastMarchingStoppingCriterionBase() {}

  OutputPixelType m_PreviousValue;
  OutputPixelType m_CurrentValue;

 private:
  FastMarchingStoppingCriterionBase( const Self& );
  void operator = ( const Self& );
 };
 }
 #endif

