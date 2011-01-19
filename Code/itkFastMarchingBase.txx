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

#ifndef __itkFastMarchingBase_txx
#define __itkFastMarchingBase_txx

#include "itkImage.h"
#include "itkImageToImageFilter.h"

#include "itkMesh.h"
#include "itkMeshToMeshFilter.h"

#include "itkPriorityQueueContainer.h"

namespace itk
{
// -----------------------------------------------------------------------------
template< class TTraits >
FastMarchingBase< TTraits >::
FastMarchingBase()
  {
  m_Heap = PriorityQueueType::New();
  m_SpeedConstant = 1.;
  m_InverseSpeed = 1.;
  m_NormalizationFactor = 1.;
  m_StoppingValue = NumericTraits< OutputPixelType >::Zero;
  m_TargetReachedValue = NumericTraits< OutputPixelType >::Zero;
  m_TopologyCheck = None;
  m_TargetCondition = NoTargets;
  m_NumberOfTargetsToBeReached = 0;
  m_LargeValue = NumericTraits< OutputPixelType >::max();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
FastMarchingBase< TTraits >::
~FastMarchingBase()
  {
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
void
FastMarchingBase< TTraits >::
PrintSelf( std::ostream & os, Indent indent ) const
  {
  Superclass::PrintSelf( os, indent );

  os << indent << "Speed constant: " << m_SpeedConstant << std::endl;
  os << indent << "Stopping value: " << m_StoppingValue << std::endl;
  os << indent << "Normalization Factor: " << m_NormalizationFactor << std::endl;
  }

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
void
FastMarchingBase< TTraits >::
SetAliveNodes( NodeContainerType iNodes )
  {
  m_AliveNodes = iNodes;
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
void
FastMarchingBase< TTraits >::
AddAliveNode( NodeType iNode, OutputPixelType iValue )
  {
  m_AliveNodes.push_back( NodePairType( iNode, iValue ) );
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
void
FastMarchingBase< TTraits >::
SetTrialNodes( NodeContainerType iNodes )
  {
  m_TrialNodes = iNodes;
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
void
FastMarchingBase< TTraits >::
AddTrialNode( NodeType iNode, OutputPixelType iValue )
  {
  m_TrialNodes.push_back( NodePairType( iNode, iValue ) );
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
void
FastMarchingBase< TTraits >::
SetForbiddenNodes( std::vector< NodeType > iNodes )
  {
  m_ForbiddenNodes = iNodes;
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
void
FastMarchingBase< TTraits >::
AddForbiddenNode( NodeType iNode )
  {
  m_ForbiddenNodes.push_back( iNode );
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
void
FastMarchingBase< TTraits >::
SetTargetNodes( std::vector< NodeType > iNodes )
  {
  m_TargetNodes = iNodes;
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
void
FastMarchingBase< TTraits >::
AddTargetNode( NodeType iNode )
  {
  m_TargetNodes.push_back( iNode );
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
bool
FastMarchingBase< TTraits >::
CheckTargetCondition( NodeType iNode )
  {
  // Only check for reached targets if the mode is not NoTargets and
  // there is at least one TargetPoint.
  if ( m_TargetCondition != NoTargets &&  !m_TargetNodes.empty() )
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
        return ( m_ReachedTargetNodes.size() == m_NumberOfTargetsToBeReached );
        }
      ++pointsIter;
      }
    }
  return false;
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
void
FastMarchingBase< TTraits >::
Initialize()
  {
  if( m_StoppingValue < 0. )
    {
    itkExceptionMacro( <<"Stopping Value is null or negative" );
    }
  if( m_TargetCondition != NoTargets )
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
        <<"Number of target nodes to be reached is above the provided number of target nodes" );
      }
    m_ReachedTargetNodes.clear();
    }
  if( m_NormalizationFactor < vnl_math::eps )
    {
    itkExceptionMacro( <<"Normalization Factor is null or negative" );
    }

  // make sure the heap is empty
  while ( !m_Heap->Empty() )
    {
    m_Heap->Pop();
    }

  InitializeOutput();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits >
void
FastMarchingBase< TTraits >::
GenerateData()
  {
  Initialize();
  try
    {
    double newProgress = 0.;
    double oldProgress = 0.;

    while( !m_Heap->Empty() )
      {
      PriorityQueueElementType element = m_Heap->Peek();
      m_Heap->Pop();

      NodeType current_node = element.m_Element;
      OutputPixelType current_value = element.m_Priority;

      // is this node already alive ?
      if( this->GetLabelValueForGivenNode( current_node ) != Alive )
        {
        this->CheckTopology( current_node );

        if( current_value > m_StoppingValue )
          {
          m_TargetReachedValue = m_StoppingValue;
          this->UpdateProgress(1.0);
          break;
          }

        // set this node as alive
        this->SetLabelValueForGivenNode( current_node, Alive );

        // update its neighbors
        this->UpdateNeighbors( current_node );

        if( !CheckTargetCondition( current_node ) )
          {
          // Send events every certain number of points.
          newProgress = static_cast< double >( current_value ) /
            static_cast< double >( m_StoppingValue );

          if ( newProgress - oldProgress > 0.01 )
            {
            this->UpdateProgress(newProgress);
            oldProgress = newProgress;
            }
          }
        else
          {
          m_TargetReachedValue = current_value;
          this->UpdateProgress(1.0);
          break;
          }
        }

      }
    }
  catch ( ProcessAborted & )
    {
    // User aborted filter excecution Here we catch an exception thrown by the
    // progress reporter and rethrow it with the correct line number and file
    // name. We also invoke AbortEvent in case some observer was interested on
    // it.
    //
    // RELEASE MEMORY!!!
    while( !m_Heap->Empty() )
      {
      m_Heap->Pop();
      }

    throw ProcessAborted(__FILE__, __LINE__);
    }

  // let's release some useless memory...
  while( !m_Heap->Empty() )
    {
    m_Heap->Pop();
    }
  }
// -----------------------------------------------------------------------------

} // end of namespace itk

#endif
