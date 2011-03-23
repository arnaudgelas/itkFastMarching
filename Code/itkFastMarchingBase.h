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

#include "itkIntTypes.h"
#include "itkFastMarchingStoppingCriterionBase.h"

//#include "itkPriorityQueueContainer.h"
#include <queue>

namespace itk
{
/**
 * \class FastMarchingBase
 * \brief Abstract class to solve an Eikonal based-equation using Fast Marching
 * Method.
 *
 * Fast marching solves an Eikonal equation where the speed is always
 * non-negative and depends on the position only. Starting from an
 * initial position on the front, fast marching systematically moves the
 * front forward one node at a time.
 *
 * Updates are preformed using an entropy satisfy scheme where only
 * "upwind" neighborhoods are used. This implementation of Fast Marching
 * uses a std::priority_queue to locate the next proper node to
 * update.
 *
 * Fast Marching sweeps through N grid points in (N log N) steps to obtain
 * the arrival time value as the front propagates through the domain.
 *
 * Implementation of this class is based on Chapter 8 of
 * "Level Set Methods and Fast Marching Methods", J.A. Sethian,
 * Cambridge Press, Second edition, 1999.
 *
 * The initial front is specified by two containers:
 * \li one containing the known nodes (Alive Nodes: nodes that are already part
 * of the object),
 * \li one containing the trial nodes (Trial Nodes: nodes that are considered
 * for inclusion).
 *
 * In order for the filter to evolve, at least some trial nodes must be
 * specified. These can for instance be specified as the layer of
 * nodes around the alive ones.
 *
 * The algorithm is terminated early by setting an appropriate stopping
 * criterion, or if there are no more nodes to process.
 *
 * \todo In the current implemenation, std::priority_queue only allows
 * taking nodes out from the front and putting nodes in from the back.
 * Use itk::PriorityQueueContainer instead.
 *
 * \sa FastMarchingStoppingCriterionBase
*/
template< class TTraits >
class FastMarchingBase : public TTraits::SuperclassType
  {
public:
  typedef TTraits                               Traits;
  typedef typename Traits::SuperclassType       SuperclassType;

  typedef FastMarchingBase            Self;
  typedef SuperclassType              Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;


  typedef typename Traits::InputDomainType        InputDomainType;
  typedef typename Traits::InputDomainPointer     InputDomainPointer;
  typedef typename Traits::InputPixelType         InputPixelType;

  typedef typename Traits::OutputDomainType       OutputDomainType;
  typedef typename Traits::OutputDomainPointer    OutputDomainPointer;
  typedef typename Traits::OutputPixelType        OutputPixelType;

  typedef typename Traits::NodeType                 NodeType;
  typedef typename Traits::NodePairType             NodePairType;
  typedef typename Traits::NodePairContainerType    NodePairContainerType;
  typedef typename Traits::NodePairContainerPointer NodePairContainerPointer;
  typedef typename Traits::NodePairContainerConstIterator
    NodePairContainerConstIterator;

  typedef typename Traits::NodeContainerType        NodeContainerType;
  typedef typename Traits::NodeContainerPointer     NodeContainerPointer;
  typedef typename Traits::NodeContainerConstIterator
    NodeContainerConstIterator;

  typedef FastMarchingStoppingCriterionBase< NodeType, OutputPixelType >
    StoppingCriterionType;
  typedef typename StoppingCriterionType::Pointer StoppingCriterionPointer;

  /*
  typedef long ElementIdentifier;

  typedef MinPriorityQueueElementWrapper< NodeType,
    OutputPixelType,
    ElementIdentifier > PriorityQueueElementType;

  typedef PriorityQueueContainer< PriorityQueueElementType,
    PriorityQueueElementType,
    OutputPixelType,
    ElementIdentifier > PriorityQueueType;
  typedef typename PriorityQueueType::Pointer PriorityQueuePointer;
  */


  /** \enum LabelType Fast Marching algorithm nodes types. */
  enum LabelType {
    /** \c Far represent far away nodes*/
    Far = 0,
    /** \c Alive represent nodes which have already been processed*/
    Alive,
    /** \c Trial represent nodes within a narrowband of the propagating front */
    Trial,
    /** \c InitialTrial represent nodes from where the propagation is initiated */
    InitialTrial,
    /** \c Forbidden represent nodes where the front can not propagate */
    Forbidden,
    /** \c Topology represent trial nodes but their inclusion would have violated
    topology checks. */
    Topology };

  /** \enum TopologyCheckType */
  enum TopologyCheckType {
    /** \c None */
    None = 0,
    /** \c NoHandles */
    NoHandles,
    /** \c Strict */
    Strict };

  /** Set/Get boolean macro indicating whether the user wants to check topology. */
  itkSetMacro( TopologyCheck, TopologyCheckType );
  itkGetConstReferenceMacro( TopologyCheck, TopologyCheckType );

  itkSetObjectMacro( TrialPoints, NodePairContainerType );
  itkGetObjectMacro( TrialPoints, NodePairContainerType );

  itkSetObjectMacro( AlivePoints, NodePairContainerType );
  itkGetObjectMacro( AlivePoints, NodePairContainerType );

  itkSetObjectMacro( ProcessedPoints, NodePairContainerType );
  itkGetObjectMacro( ProcessedPoints, NodePairContainerType );

  itkSetObjectMacro( ForbiddenPoints, NodeContainerType );
  itkGetObjectMacro( ForbiddenPoints, NodeContainerType );

  /** \brief Set/Get the Stopping Criterion */
  itkGetObjectMacro( StoppingCriterion, StoppingCriterionType );
  itkSetObjectMacro( StoppingCriterion, StoppingCriterionType );

  itkGetMacro( SpeedConstant, double );
  itkSetMacro( SpeedConstant, double );

  itkGetMacro( TargetReachedValue, OutputPixelType );

  /** Set the Collect Points flag. Instrument the algorithm to collect
  * a container of all nodes which it has visited. Useful for
  * creating Narrowbands for level set algorithms that supports
  * narrow banding. */
  itkSetMacro(CollectPoints, bool);

  /** Get thConste Collect Points flag. */
  itkGetConstReferenceMacro(CollectPoints, bool);
  itkBooleanMacro(CollectPoints);


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

  NodePairContainerPointer  m_TrialPoints;
  NodePairContainerPointer  m_AlivePoints;
  NodePairContainerPointer  m_ProcessedPoints;
  NodeContainerPointer      m_ForbiddenPoints;

  bool m_CollectPoints;

  //PriorityQueuePointer m_Heap;
  typedef std::vector< NodePairType >   HeapContainerType;
  typedef std::greater< NodePairType >  NodeComparerType;

  typedef std::priority_queue< NodeType, HeapContainerType, NodeComparerType >
    PriorityQueueType;

  PriorityQueueType m_Heap;

  TopologyCheckType m_TopologyCheck;

  virtual IdentifierType GetTotalNumberOfNodes() const = 0;

  virtual const OutputPixelType GetOutputValue( OutputDomainType* oDomain,
                                         const NodeType& iNode ) const = 0;

  virtual void SetOutputValue( OutputDomainType* oDomain,
                              const NodeType& iNode,
                              const OutputPixelType& iValue ) = 0;

  /** \brief Get the LabelType Value for a given node
    \param[in] iNode
    \return its label value  */
  virtual const unsigned char
  GetLabelValueForGivenNode( const NodeType& iNode ) const = 0;

  /** \brief Set the Label Value for a given node
    \param[in] iNode
    \param[in] iLabel */
  virtual void SetLabelValueForGivenNode( const NodeType& iNode,
                                         const LabelType& iLabel ) = 0;

  /** \brief Update neighbors to a given node
    \param[in] oDomain
    \param[in] iNode
  */
  virtual void UpdateNeighbors( OutputDomainType* oDomain,
                               const NodeType& iNode ) = 0;

  /** \brief Update value for a given node
    \param[in] oDomain
    \param[in] iNode
    */
  virtual void UpdateValue( OutputDomainType* oDomain,
                           const NodeType& iNode ) = 0;

  /** \brief Check if the current node violate topological criterion.
    \param[in] oDomain
    \param[in] iNode
   */
  virtual bool CheckTopology( OutputDomainType* oDomain,
                             const NodeType& iNode ) = 0;

  /** \brief   */
  void Initialize( OutputDomainType* oDomain );

  /**    */
  virtual void InitializeOutput( OutputDomainType* oDomain ) = 0;

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
