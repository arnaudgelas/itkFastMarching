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
template< class TTraits, class TCriterion >
FastMarchingBase< TTraits, TCriterion >::
FastMarchingBase()
  {
  this->ProcessObject::SetNumberOfRequiredInputs(0);

  //m_Heap = PriorityQueueType::New();
  m_SpeedConstant = 1.;
  m_InverseSpeed = -1.;
  m_NormalizationFactor = 1.;
  m_TargetReachedValue = NumericTraits< OutputPixelType >::Zero;
  m_TopologyCheck = None;
  m_LargeValue = NumericTraits< OutputPixelType >::max();
  m_TopologyValue = -1. * std::numeric_limits< OutputPixelType >::epsilon();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
FastMarchingBase< TTraits, TCriterion >::
~FastMarchingBase()
  {
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
void
FastMarchingBase< TTraits, TCriterion >::
PrintSelf( std::ostream & os, Indent indent ) const
  {
  Superclass::PrintSelf( os, indent );

  os << indent << "Speed constant: " << m_SpeedConstant << std::endl;
  os << indent << "Normalization Factor: " << m_NormalizationFactor << std::endl;
  }

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
void
FastMarchingBase< TTraits, TCriterion >::
SetAliveNodes( NodeContainerType iNodes )
  {
  m_AliveNodes = iNodes;
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
void
FastMarchingBase< TTraits, TCriterion >::
AddAliveNode( const NodeType& iNode, const OutputPixelType& iValue )
  {
  m_AliveNodes.push_back( NodePairType( iNode, iValue ) );
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
void
FastMarchingBase< TTraits, TCriterion >::
AddAliveNode( const NodePairType& iPair )
  {
  m_AliveNodes.push_back( iPair );
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
void
FastMarchingBase< TTraits, TCriterion >::
SetTrialNodes( NodeContainerType iNodes )
  {
  m_TrialNodes = iNodes;
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
void
FastMarchingBase< TTraits, TCriterion >::
AddTrialNode( const NodeType& iNode, const OutputPixelType& iValue )
  {
  m_TrialNodes.push_back( NodePairType( iNode, iValue ) );
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
void
FastMarchingBase< TTraits, TCriterion >::
AddTrialNode( const NodePairType& iPair )
  {
  m_TrialNodes.push_back( iPair );
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
void
FastMarchingBase< TTraits, TCriterion >::
SetForbiddenNodes( std::vector< NodeType > iNodes )
  {
  m_ForbiddenNodes = iNodes;
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
void
FastMarchingBase< TTraits, TCriterion >::
AddForbiddenNode( const NodeType& iNode )
  {
  m_ForbiddenNodes.push_back( iNode );
  this->Modified();
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
void
FastMarchingBase< TTraits, TCriterion >::
Initialize( OutputDomainType* oDomain )
  {
  if( m_StoppingCriterion.IsNull() )
    {
    itkExceptionMacro( <<"No Stopping Criterion Set" );
    }
  if( m_NormalizationFactor < vnl_math::eps )
    {
    itkExceptionMacro( <<"Normalization Factor is null or negative" );
    }
  if( m_SpeedConstant < vnl_math::eps )
    {
    itkExceptionMacro( <<"SpeedConstant is null or negative" );
    }

  // make sure the heap is empty
  while ( !m_Heap.empty() )
    {
    m_Heap.pop();
    }
  /*
  while ( !m_Heap->Empty() )
    {
    m_Heap->Pop();
    }
  */

  InitializeOutput( oDomain );
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< class TTraits, class TCriterion >
void
FastMarchingBase< TTraits, TCriterion >::
GenerateData()
  {
  OutputDomainType* output = this->GetOutput();

  Initialize( output );

  OutputPixelType current_value = 0.;

  this->UpdateProgress(0.0);   // Send first progress event

  try
    {
    //double newProgress = 0.;
    //double oldProgress = 0.;

    //while( !m_Heap->Empty() )
    while( !m_Heap.empty() )
      {
      //PriorityQueueElementType element = m_Heap->Peek();
      //m_Heap->Pop();
      //
      //NodeType current_node = element.m_Element;
      //OutputPixelType current_value = element.m_Priority;


      NodePairType current_node_pair = m_Heap.top();
      m_Heap.pop();

      NodeType current_node = current_node_pair.first;
      current_value = this->GetOutputValue( output, current_node );

      if( current_value == current_node_pair.second )
        {
        // is this node already alive ?
        if( this->GetLabelValueForGivenNode( current_node ) != Alive )
          {
          m_StoppingCriterion->SetCurrentNode( current_node );
          m_StoppingCriterion->SetCurrentValue( current_value );

          if( m_StoppingCriterion->IsSatisfied() )
            {
            break;
            }

          if( this->CheckTopology( output, current_node ) )
            {
            // set this node as alive
            this->SetLabelValueForGivenNode( current_node, Alive );

            // update its neighbors
            this->UpdateNeighbors( output, current_node );


            // Send events every certain number of points.
            /*
            newProgress = static_cast< double >( current_value ) /
              static_cast< double >( m_StoppingValue );

            if ( newProgress - oldProgress > 0.01 )
              {
              this->UpdateProgress(newProgress);
              oldProgress = newProgress;
              }*/
            }
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
    while( !m_Heap.empty() )
      {
      m_Heap.pop();
      }
    /*while( !m_Heap->Empty() )
      {
      m_Heap->Pop();
      }*/

    throw ProcessAborted(__FILE__, __LINE__);
    }

  m_TargetReachedValue = current_value;
  this->UpdateProgress(1.0);

  // let's release some useless memory...
  while( !m_Heap.empty() )
    {
    m_Heap.pop();
    }
  /*while( !m_Heap->Empty() )
    {
    m_Heap->Pop();
    }*/
  }
// -----------------------------------------------------------------------------

} // end of namespace itk

#endif
