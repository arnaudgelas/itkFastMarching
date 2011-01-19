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

#ifndef __itkFastMarchingBase_h
#define __itkFastMarchingBase_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"

#include "itkMesh.h"
#include "itkMeshToMeshFilter.h"

#include "itkPriorityQueueContainer.h"

namespace itk
{
template< class TInputDomain,
          class TNode,
          class TOutputDomain,
          class TSuperclass >
class FastMarchingTraits
  {
public:
  typedef TInputDomain InputDomainType;
  typedef typename InputDomainType::Pointer InputDomainPointer;
  typedef typename InputDomainType::PixelType InputPixelType;

  typedef TNode NodeType;
  typedef TOutputDomain OutputDomainType;
  typedef typename OutputDomainType::Pointer OutputDomainPointer;
  typedef typename OutputDomainType::PixelType OutputPixelType;
  typedef TSuperclass SuperclassType;

  // Here concept checking: make sure TValue is scalar!
  };

template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
class ImageFastMarchingTraits :
    public FastMarchingTraits<
      Image< TInputPixel, VDimension >,
      Index< VDimension >,
      Image< TOutputPixel, VDimension >,
      ImageToImageFilter< Image< TInputPixel, VDimension >,
                          Image< TOutputPixel, VDimension > >
    >
  {

  };

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
class MeshFastMarchingTraits :
    public FastMarchingTraits<
      Mesh< TInputPixel, VDimension, TInputMeshTraits >,
      typename TInputMeshTraits::PointIdentifier,
      Mesh< TOutputPixel, VDimension, TOutputMeshTraits >,
      MeshToMeshFilter< Mesh< TInputPixel, VDimension, TInputMeshTraits >,
                        Mesh< TOutputPixel, VDimension, TOutputMeshTraits > > >
  {};


/**
 * \class FastMarchingBase
 * \brief Solve an Eikonal equation (see equation below) using Fast Marching
 *
 * Fast marching solves an Eikonal equation where the speed is always
 * non-negative and depends on the position only. Starting from an
 * initial position on the front, fast marching systematically moves the
 * front forward one grid point at a time.
 *
 * Updates are preformed using an entropy satisfy scheme where only
 * "upwind" neighborhoods are used.
 *
 * Fast Marching sweeps through N grid points in (N log N) steps to obtain
 * the arrival time value as the front propagates through the grid.
 *
 * Implementation of this class is based on Chapter 8 of
 * "Level Set Methods and Fast Marching Methods", J.A. Sethian,
 * Cambridge Press, Second edition, 1999.
 *
*/
template< class TTraits >
class FastMarchingBase : public TTraits::SuperclassType
  {
public:
  typedef TTraits                               Traits;
  typedef typename Traits::SuperclassType     SuperclassType;

  typedef FastMarchingBase            Self;
  typedef SuperclassType              Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;


  typedef typename Traits::InputDomainType     InputDomainType;
  typedef typename Traits::InputDomainPointer  InputDomainPointer;
  typedef typename Traits::InputPixelType      InputPixelType;

  typedef typename Traits::OutputDomainType     OutputDomainType;
  typedef typename Traits::OutputDomainPointer  OutputDomainPointer;
  typedef typename Traits::OutputPixelType      OutputPixelType;

  typedef typename Traits::NodeType NodeType;
  typedef long ElementIdentifier;

  typedef MinPriorityQueueElementWrapper< NodeType,
    OutputPixelType,
    ElementIdentifier > PriorityQueueElementType;

  typedef PriorityQueueContainer< PriorityQueueElementType,
    PriorityQueueElementType,
    OutputPixelType,
    ElementIdentifier > PriorityQueueType;
  typedef typename PriorityQueueType::Pointer PriorityQueuePointer;

  /** Enum of Fast Marching algorithm point types.
   * Far represent far away points;
   * Trial represent points within a narrowband of the propagating front;
   * Alive represent points which have already been processed.
   * Forbidden represent points where the front can not propagate
   */
  enum LabelType { Far = 0,
                   Alive,
                   Trial,
                   InitialTrial,
                   Forbidden };

  enum TopologyCheckType { None = 0,
                           NoHandles,
                           Strict };

  enum TargetConditionType { NoTargets = 0,
                             OneTarget,
                             SomeTargets,
                             AllTargets };

  /** Set/Get boolean macro indicating whether the user wants to check topology. */
  itkSetMacro( TopologyCheck, TopologyCheckType );
  itkGetConstReferenceMacro( TopologyCheck, TopologyCheckType );

  /** Set/Get boolean macro indicating whether the user wants to check topology. */
  itkSetMacro( TargetCondition, TargetConditionType );
  itkGetConstReferenceMacro( TargetCondition, TargetConditionType );

  typedef std::pair< NodeType, OutputPixelType > NodePairType;
  typedef std::vector< std::pair< NodeType, OutputPixelType > > NodeContainerType;
  typedef typename NodeContainerType::iterator NodeContainerIterator;

  virtual void SetAliveNodes( NodeContainerType iNodes )
    {
    m_AliveNodes = iNodes;
    this->Modified();
    }

  virtual void AddAliveNode( NodeType iNode, OutputPixelType iValue )
    {
    m_AliveNodes.push_back( NodePairType( iNode, iValue ) );
    this->Modified();
    }

  virtual void SetTrialNodes( NodeContainerType iNodes )
    {
    m_TrialNodes = iNodes;
    this->Modified();
    }

  virtual void AddTrialNode( NodeType iNode, OutputPixelType iValue )
    {
    m_TrialNodes.push_back( NodePairType( iNode, iValue ) );
    this->Modified();
    }

  virtual void SetForbiddenNodes( std::vector< NodeType > iNodes )
    {
    m_ForbiddenNodes = iNodes;
    this->Modified();
    }

  virtual void AddForbiddenNode( NodeType iNode )
    {
    m_ForbiddenNodes.push_back( iNode );
    this->Modified();
    }

  virtual void SetTargetNodes( std::vector< NodeType > iNodes )
    {
    m_TargetNodes = iNodes;
    this->Modified();
    }

  virtual void AddTargetNode( NodeType iNode )
    {
    m_TargetNodes.push_back( iNode );
    this->Modified();
    }


protected:

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
    m_NumberOfTargetsToBeReached = 1;
    }
  virtual ~FastMarchingBase() {}

  double m_SpeedConstant;
  double m_InverseSpeed;
  double m_NormalizationFactor;
  OutputPixelType m_StoppingValue;
  OutputPixelType m_TargetReachedValue;

  NodeContainerType m_TrialNodes;
  NodeContainerType m_AliveNodes;
  std::vector< NodeType > m_ForbiddenNodes;
  std::vector< NodeType > m_TargetNodes;
  std::vector< NodeType > m_ReachedTargetNodes;
  size_t m_NumberOfTargetsToBeReached;

  PriorityQueuePointer m_Heap;

  TopologyCheckType m_TopologyCheck;
  TargetConditionType m_TargetCondition;

  virtual LabelType GetLabelValueForGivenNode( NodeType iNode ) = 0;
  virtual void SetLabelValueForGivenNode( NodeType iNode, LabelType iLabel ) = 0;
  virtual void UpdateNeighbors( NodeType iNode ) = 0;
  virtual void CheckTopology( NodeType iNode ) = 0;

  bool CheckTargetCondition( NodeType iNode )
    {
    // Only check for reached targets if the mode is not NoTargets and
    // there is at least one TargetPoint.
    if ( m_TargetCondition != NoTargets &&  !m_TargetNodes.empty() )
      {
      typename std::vector< NodeType >::const_iterator
          pointsIter = m_TargetNodes.begin();
      typename std::vector< NodeType >::const_iterator
          pointsEnd = m_TargetNodes.end();

      switch( m_TargetCondition )
        {
        default:
        case OneTarget:
          {
          while( pointsIter != pointsEnd )
            {
            if ( *pointsIter == iNode )
              {
              return true;
              }
            ++pointsIter;
            }
          break;
          }
        case SomeTargets:
          {
          while( pointsIter != pointsEnd )
            {
            if ( *pointsIter == iNode )
              {
              this->m_ReachedTargetNodes.push_back( iNode );
              return ( m_ReachedTargetNodes.size() == m_NumberOfTargetsToBeReached );
              }
            ++pointsIter;
            }
          break;
          }
        case AllTargets:
          {
          while( pointsIter != pointsEnd )
            {
            if ( *pointsIter == iNode )
              {
              this->m_ReachedTargetNodes.push_back( iNode );

              return( m_ReachedTargetNodes.size() == m_TargetNodes.size() );
              }
            ++pointsIter;
            }
          break;
          }
        }
      }
    return false;
    }


  virtual void GenerateData()
    {
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
      ProcessAborted e(__FILE__, __LINE__);
      e.SetDescription("Process aborted.");
      e.SetLocation(ITK_LOCATION);
      throw e;
      }
    }

  void PrintSelf(std::ostream & os, Indent indent) const
    {
    Superclass::PrintSelf( os, indent );
    }



private:
  FastMarchingBase( const Self& );
  void operator = ( const Self& );
  };


}

#endif
