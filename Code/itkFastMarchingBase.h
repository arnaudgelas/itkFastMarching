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
                   Outside };

  enum TopologyCheckType { None = 0,
                           NoHandles,
                           Strict };

  /** Set/Get boolean macro indicating whether the user wants to check topology. */
  itkSetMacro( TopologyCheck, TopologyCheckType );
  itkGetConstReferenceMacro( TopologyCheck, TopologyCheckType );

  typedef std::vector< std::pair< NodeType, OutputPixelType > > NodeContainerType;
  typedef typename NodeContainerType::iterator NodeContainerIterator;

  virtual void SetAliveNodes( NodeContainerType iNodes ) = 0;
  virtual void AddAliveNode( NodeType iNode, OutputPixelType iValue ) = 0;

  virtual void SetTrialNodes( NodeContainerType iNodes ) = 0;
  virtual void AddTrialNode( NodeType iNode, OutputPixelType iValue ) = 0;

  virtual void SetForbiddenNodes( NodeType iNodes ) = 0;
  virtual void AddForbiddenNode( NodeType iNode ) = 0;


protected:

  FastMarchingBase() {}
  virtual ~FastMarchingBase() {}

  double m_SpeedConstand;
  double m_InverseSpeed;
  double m_NormalizationFactor;
  OutputPixelType m_StoppingValue;

  PriorityQueuePointer m_Heap;

  TopologyCheckType m_TopologyCheck;

  virtual LabelType GetLabelValueForGivenNode( NodeType iNode ) = 0;
  virtual void SetLabelValueForGivenNode( NodeType iNode, LabelType iLabel ) = 0;
  virtual void UpdateNeighbors( NodeType iNode ) = 0;
  virtual void CheckTopology( NodeType iNode ) = 0;

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
            this->UpdateProgress(1.0);
            break;
            }

          // set this node as alive
          this->SetLabelValueForGivenNode( current_node, Alive );

          // update its neighbors
          this->UpdateNeighbors( current_node );

          // Send events every certain number of points.
          newProgress = static_cast< double >( current_value ) /
            static_cast< double >( m_StoppingValue );
          if ( current_value - oldProgress > 0.01 )
            {
            this->UpdateProgress(newProgress);
            oldProgress = newProgress;
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



private:
  FastMarchingBase( const Self& );
  void operator = ( const Self& );
  };


}

#endif // ITKFASTMARCHING_H
