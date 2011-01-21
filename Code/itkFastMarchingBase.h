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

#include "itkFastMarchingStoppingCriterionBase.h"

#include "itkPriorityQueueContainer.h"

namespace itk
{
/**  \class FastMarchingTraits
  \brief
  */
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

/**  \class ImageFastMarchingTraits
  \brief
*/
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

/**  \class MeshFastMarchingTraits
  \brief
*/
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

  typedef FastMarchingStoppingCriterionBase< Self > StoppingCriterionType;
  typedef typename StoppingCriterionType::Pointer   StoppingCriterionPointer;

  /** Enum of Fast Marching algorithm point types.
   * Far represent far away points;
   * Trial represent points within a narrowband of the propagating front;
   * Alive represent points which have already been processed.
   * Forbidden represent points where the front can not propagate.
   * Topology points were trial points but their inclusion
   * would have violated topology checks.
   */
  enum LabelType { Far = 0,
                   Alive,
                   Trial,
                   InitialTrial,
                   Forbidden,
                   Topology };

  enum TopologyCheckType { None = 0,
                           NoHandles,
                           Strict };

  /** Set/Get boolean macro indicating whether the user wants to check topology. */
  itkSetMacro( TopologyCheck, TopologyCheckType );
  itkGetConstReferenceMacro( TopologyCheck, TopologyCheckType );

  typedef std::pair< NodeType, OutputPixelType > NodePairType;
  typedef std::vector< std::pair< NodeType, OutputPixelType > > NodeContainerType;
  typedef typename NodeContainerType::const_iterator NodeContainerConstIterator;
  typedef typename NodeContainerType::iterator NodeContainerIterator;

  /** \brief Set Alive nodes representing the initial front. */
  virtual void SetAliveNodes( NodeContainerType iNodes );

  /** \brief Add one node as alive to the initial front.
    \param[in] iNode
    \param[in] ivalue */
  virtual void AddAliveNode( NodeType iNode, OutputPixelType iValue );

  /** \brief Set Trial Nodes */
  virtual void SetTrialNodes( NodeContainerType iNodes );
  virtual void AddTrialNode( NodeType iNode, OutputPixelType iValue );

  /** \brief Set Forbidden Nodes */
  virtual void SetForbiddenNodes( std::vector< NodeType > iNodes );
  virtual void AddForbiddenNode( NodeType iNode );


protected:

  /** \brief Constructor */
  FastMarchingBase();

  /** \brief Destructor */
  virtual ~FastMarchingBase();

  StoppingCriterionPointer m_StoppingCriterion;

  double m_SpeedConstant;
  double m_InverseSpeed;
  double m_NormalizationFactor;
  OutputPixelType m_TargetReachedValue;
  OutputPixelType m_LargeValue;
  OutputPixelType m_TopologyValue;

  NodeContainerType m_TrialNodes;
  NodeContainerType m_AliveNodes;
  std::vector< NodeType > m_ForbiddenNodes;

  PriorityQueuePointer m_Heap;

  TopologyCheckType m_TopologyCheck;

  /** \brief Get the Label Value for a given node
    \param[in] iNode
    \return its label value  */
  virtual char GetLabelValueForGivenNode( NodeType iNode ) = 0;

  /** \brief Set the Label Value for a given node
    \param[in] iNode
    \param[in] iLabel */
  virtual void SetLabelValueForGivenNode( NodeType iNode, LabelType iLabel ) = 0;

  /** \brief Update neighbors to a given node
    \param[in] iNode
  */
  virtual void UpdateNeighbors( NodeType iNode ) = 0;

  /** \brief Update value for a given node
    \param[in] iNode
    */
  virtual void UpdateValue( NodeType iNode ) = 0;

  /**    */
  virtual bool CheckTopology( NodeType iNode ) = 0;

  /**    */
  void Initialize();

  /**    */
  virtual void InitializeOutput() = 0;

  /**    */
  void GenerateData();

  /**    */
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  FastMarchingBase( const Self& );
  void operator = ( const Self& );
  };

}

#include "itkFastMarchingBase.txx"
#endif
