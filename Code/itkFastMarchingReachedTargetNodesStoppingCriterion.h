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
  template< class TNode,
           typename TValue >
  class FastMarchingReachedTargetNodesStoppingCriterion :
      public FastMarchingStoppingCriterionBase< TNode, TValue >
    {
  public:
    typedef FastMarchingReachedTargetNodesStoppingCriterion Self;
    typedef FastMarchingStoppingCriterionBase< TNode, TValue > Superclass;
    typedef SmartPointer< Self >              Pointer;
    typedef SmartPointer< const Self >        ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(FastMarchingReachedTargetNodesStoppingCriterion,
                 FastMarchingStoppingCriterionBase );

    typedef typename Superclass::ValueType  ValueType;
    typedef typename Superclass::NodeType   NodeType;

    enum TargetConditionType { OneTarget = 1,
                               SomeTargets,
                               AllTargets };

    /** Set/Get boolean macro indicating whether the user wants to check topology. */
    void SetTargetCondition( TargetConditionType iCondition )
      {
      m_TargetCondition = iCondition;
      m_Initialized = false;
      this->Modified();
      }

    itkGetConstReferenceMacro( TargetCondition, TargetConditionType );

    itkSetMacro( TargetOffset, ValueType );
    itkGetMacro( TargetOffset, ValueType );

    void SetNumberOfTargetsToBeReached( size_t iN )
      {
      m_NumberOfTargetsToBeReached = iN;
      m_Initialized = false;
      this->Modified();
      }

    /** \brief Set Target Nodes*/
    virtual void SetTargetNodes( const std::vector< NodeType >& iNodes )
      {
      m_TargetNodes = iNodes;
      m_Initialized = false;
      this->Modified();
      }

    void SetCurrentNode( const NodeType& iNode )
      {
      if( !m_Initialized )
        {
        Initialize();
        }

      if( !m_Satisfied )
        {
        // Only check for reached targets if the mode is not NoTargets and
        // there is at least one TargetPoint.
        if ( !m_TargetNodes.empty() )
          {
          typename std::vector< NodeType >::const_iterator
              pointsIter = m_TargetNodes.begin();
          typename std::vector< NodeType >::const_iterator
              pointsEnd = m_TargetNodes.end();

          while( pointsIter != pointsEnd )
            {
            if ( *pointsIter == iNode )
              {
              this->m_ReachedTargetNodes.push_back( iNode );
              m_Satisfied =
                  ( m_ReachedTargetNodes.size() == m_NumberOfTargetsToBeReached );
              break;
              }
            ++pointsIter;
            }
          if( m_Satisfied )
            {
            m_StoppingValue = this->m_CurrentValue + m_TargetOffset;
            }
          }
        else
          {
          m_Satisfied = false;
          }
        }
      }

    bool IsSatisfied() const
      {
      return m_Satisfied && ( this->m_CurrentValue >= m_StoppingValue );
      }

    const std::string GetDescription() const
      {
      return "Target Nodes Reached with possible overshoot";
      }

  protected:
    FastMarchingReachedTargetNodesStoppingCriterion() : Superclass()
      {
      m_TargetCondition = AllTargets;
      m_TargetOffset = NumericTraits< ValueType >::Zero;
      m_StoppingValue = NumericTraits< ValueType >::Zero;
      m_Satisfied = false;
      m_Initialized = false;
      }
    ~FastMarchingReachedTargetNodesStoppingCriterion() {}

    TargetConditionType m_TargetCondition;
    std::vector< NodeType > m_TargetNodes;
    std::vector< NodeType > m_ReachedTargetNodes;
    size_t m_NumberOfTargetsToBeReached;
    ValueType m_TargetOffset;
    ValueType m_StoppingValue;
    bool m_Satisfied;
    bool m_Initialized;

    void Initialize()
      {
      if( m_TargetCondition == OneTarget )
        {
        m_NumberOfTargetsToBeReached = 1;
        }
      if( m_TargetCondition == AllTargets )
        {
        m_NumberOfTargetsToBeReached = m_TargetNodes.size();
        }
      if( m_NumberOfTargetsToBeReached < 1 )
        {
        itkExceptionMacro(
              <<"Number of target nodes to be reached is null" );
        }
      if( m_NumberOfTargetsToBeReached > m_TargetNodes.size() )
        {
        itkExceptionMacro(
          <<"Number of target nodes to be reached is above the provided number of           target nodes" );
        }
      m_ReachedTargetNodes.clear();

      m_Satisfied = false;
      m_Initialized = true;
      }

  private:
    FastMarchingReachedTargetNodesStoppingCriterion( const Self& );
    void operator = ( const Self& );
    };

  template< class TImage >
  class FastMarchingImageReachedTargetNodesStoppingCriterion :
    public FastMarchingReachedTargetNodesStoppingCriterion<
      typename TImage::IndexType,
      typename TImage::PixelType >
    {
  public:
    typedef TImage ImageType;
    typedef typename ImageType::IndexType NodeType;
    typedef typename ImageType::PixelType ValueType;

    typedef FastMarchingImageReachedTargetNodesStoppingCriterion Self;
    typedef FastMarchingReachedTargetNodesStoppingCriterion< NodeType,
      ValueType > Superclass;
    typedef SmartPointer< Self >              Pointer;
    typedef SmartPointer< const Self >        ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(FastMarchingImageReachedTargetNodesStoppingCriterion,
                 FastMarchingReachedTargetNodesStoppingCriterion );

  protected:
    FastMarchingImageReachedTargetNodesStoppingCriterion() : Superclass() {}
    ~FastMarchingImageReachedTargetNodesStoppingCriterion() {}
    };
}
#endif // __itkFastMarchingThresholdStoppingCriterion_h
