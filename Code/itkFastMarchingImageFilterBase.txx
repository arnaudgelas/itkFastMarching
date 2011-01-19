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

#ifndef __itkFastMarchingImageFilterBase_txx
#define __itkFastMarchingImageFilterBase_txx

#include "itkImageRegionIterator.h"

namespace itk
{
// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
class
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
    InternalNodeStructure
  {
public:
    InternalNodeStructure( ) : m_Value( this->m_LargeValue ) {}

    NodeType        m_Node;
    OutputPixelType m_Value;
    unsigned int    m_Axis;

    bool operator< ( InternalNodeStructure& iRight )
      {
      return m_Value < iRight.m_Value;
      }
  };


// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
FastMarchingImageFilterBase()
  {
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
~FastMarchingImageFilterBase()
  {
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
typename
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
LabelType
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
GetLabelValueForGivenNode( NodeType iNode )
  {
  return m_LabelImage->GetPixel( iNode );
  }
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
void
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
SetLabelValueForGivenNode( NodeType iNode, LabelType iLabel )
  {
  m_LabelImage->SetPixel( iNode, iLabel );
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
void
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
UpdateNeighbors( NodeType iNode )
  {
  NodeType neighIndex = iNode;
  unsigned char label;

  for ( unsigned int j = 0; j < ImageDimension; j++ )
    {
    for( int s = -1; s < 2; s+= 2 )
      {
      if ( ( iNode[j] > m_StartIndex[j] ) && ( iNode[j] < m_LastIndex[j] ) )
        {
        neighIndex[j] += s;
        }
     label = m_LabelImage->GetPixel(neighIndex);

      if ( ( label != Superclass::Alive ) &&
           ( label != Superclass::InitialTrial ) &&
           ( label != Superclass::Forbidden ) )
        {
        this->UpdateValue( neighIndex );
        }
      }

    //reset neighIndex
    neighIndex[j] = iNode[j];
    }
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
void
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
UpdateValue( NodeType iNode )
  {
  NodeType neighbor_node = iNode;

  OutputPixelType neighValue;

  std::vector< InternalNodeStructure > NodesUsed( ImageDimension );

  // just to make sure the index is initialized (really cautious)
  InternalNodeStructure temp_node;
  temp_node.m_Node = iNode;

  for ( unsigned int j = 0; j < ImageDimension; j++ )
    {
    temp_node.m_Value = this->m_LargeValue;

    // find smallest valued neighbor in this dimension
    for ( int s = -1; s < 2; s = s + 2 )
      {
      neighbor_node[j] = iNode[j] + s;

      // make sure neighIndex is not outside from the image
      if ( ( neighbor_node[j] > m_LastIndex[j] ) ||
           ( neighbor_node[j] < m_StartIndex[j] ) )
        {
        continue;
        }

      if ( m_LabelImage->GetPixel( neighbor_node ) == Superclass::Alive )
        {
        neighValue =
            static_cast< OutputPixelType >( this->GetOutput()->GetPixel( neighbor_node ) );

        // let's find the minimum value given a direction j
        if ( temp_node.m_Value > neighValue )
          {
          temp_node.m_Value = neighValue;
          temp_node.m_Node = neighbor_node;
          }
        }
      } // end for ( int s = -1; s < 2; s = s + 2 )

    // put the minimum neighbor onto the heap
    temp_node.m_Axis = j;
    NodesUsed[j] = temp_node;

    // reset neighIndex
    neighbor_node[j] = iNode[j];

    } // end for ( unsigned int j = 0; j < SetDimension; j++ )

  double solution = Solve( NodesUsed );

  if ( solution < this->m_LargeValue )
    {
    // write solution to m_OutputLevelSet
    OutputPixelType outputPixel = static_cast< OutputPixelType >( solution );
    this->GetOutput()->SetPixel(index, outputPixel);

    // insert point into trial heap
    m_LabelImage->SetPixel( index, Superclass::Trial );
    //node.SetValue( outputPixel );
    //node.SetIndex( index );
    //m_TrialHeap.push(node);
    }
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
double
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
Solve( std::vector< InternalNodeStructure > iNeighbors )
{
  // sort the local list
  std::sort( iNeighbors.begin(), iNeighbors.end() );

  double oSolution = NumericTraits< double >::max();

  double aa( 0.0 );
  double bb( 0.0 );
  double cc( this->m_InverseSpeed );

  if ( this->GetInput() )
    {
    cc =
      static_cast< double >( this->GetInput()->GetPixel(index) ) /
        this->m_NormalizationFactor;
    cc = -1.0 * vnl_math_sqr(1.0 / cc);
    }

  OutputSpacingType spacing = this->GetOutput()->GetSpacing();

  double discrim = 0.;
  double value = 0.;
  double spaceFactor = 0.;
  unsigned int axis = 0;

  typename std::vector< InternalNodeStructure >::iterator
      n_it = iNeighbors.begin();

  while( n_it != iNeighbors.end() )
    {
    value = static_cast< double >( n_it->m_Value );

    if ( oSolution >= value )
      {
      axis = n_it->m_Axis;

      // spaceFactor = \frac{1}{spacing[axis]^2}
      spaceFactor = vnl_math_sqr(1.0 / spacing[axis]);

      aa += spaceFactor;
      bb += value * spaceFactor;
      cc += vnl_math_sqr(value) * spaceFactor;

      discrim = vnl_math_sqr(bb) - aa * cc;
      if ( discrim < vnl_math::eps )
        {
        // Discriminant of quadratic eqn. is negative
        itkExceptionMacro(
          <<"Discriminant of quadratic equation is negative" );
        }

      oSolution = ( vcl_sqrt(discrim) + bb ) / aa;
      }
    else
      {
      break;
      }
    ++n_it;
    }

  return oSolution;
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
void
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
CheckTopology( NodeType iNode )
  {
  (void) iNode;

  if( this->m_TopologyCheck != Superclass::None )
    {
    itkWarningMacro( << "CheckTopology has not be implemented for Dimension != 2 and != 3."
                    << "m_TopologyCheck should be set to None." );
    }
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
void FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
InitializeOutput()
  {
  OutputImageType* output = this->GetOutput();

  // allocate memory for the output buffer
  output->SetBufferedRegion( output->GetRequestedRegion() );
  output->Allocate();

  // cache some buffered region information
  //m_BufferedRegion = output->GetBufferedRegion();
  //m_StartIndex = m_BufferedRegion.GetIndex();
  //m_LastIndex = m_StartIndex + m_BufferedRegion.GetSize();

  typename OutputImageType::OffsetType offset;
  offset.Fill(1);
  m_LastIndex -= offset;

  // allocate memory for the PointTypeImage
  m_LabelImage->CopyInformation(output);
  m_LabelImage->SetBufferedRegion( output->GetBufferedRegion() );
  m_LabelImage->Allocate();

  // set all output value to infinity
  typedef ImageRegionIterator< OutputImageType > OutputIterator;

  OutputIterator outIt ( output, output->GetBufferedRegion() );

  // set all points type to FarPoint
  typedef ImageRegionIterator< LabelImageType > LabelIterator;

  LabelIterator typeIt( m_LabelImage,
                        m_LabelImage->GetBufferedRegion() );


  OutputPixelType outputPixel = this->m_LargeValue;

  outIt.GoToBegin();
  typeIt.GoToBegin();
  while( !outIt.IsAtEnd() )
    {
    outIt.Set(outputPixel);
    typeIt.Set( Superclass::Far );

    ++outIt;
    ++typeIt;
    }


  NodeType idx;

  if ( !this->m_AliveNodes.empty() )
    {
    NodeContainerConstIterator pointsIter = this->m_AliveNodes.begin();
    NodeContainerConstIterator pointsEnd = this->m_AliveNodes.end();

    while( pointsIter != pointsEnd )
      {
      // get node from alive points container
      idx = pointsIter->first;

      // check if node index is within the output level set
      //if ( m_BufferedRegion.IsInside( idx ) )
        {
        // make this an alive point
        m_LabelImage->SetPixel(idx, Superclass::Alive );

        outputPixel = pointsIter->second;
        output->SetPixel(idx, outputPixel);
        }

      ++pointsIter;
      }
    }

  if( this->m_ForbiddenNodes.empty() )
    {
    typename std::vector< NodeType >::const_iterator
        p_it = this->m_ForbiddenNodes.begin();
    typename std::vector< NodeType >::const_iterator
        p_end = this->m_ForbiddenNodes.end();

    OutputPixelType zero = NumericTraits< OutputPixelType >::Zero;

    while( p_it != p_end )
      {
      idx = *p_it;

      // check if node index is within the output level set
      //if ( m_BufferedRegion.IsInside( idx ) )
        {
        // make this an alive point
        m_LabelImage->SetPixel(idx, Superclass::Forbidden );
        output->SetPixel (idx, zero );
        }

      ++p_it;
      }
    }

  // process the input trial points
  if ( this->m_TrialNodes.empty() )
    {
    NodeContainerConstIterator pointsIter = this->m_TrialNodes.begin();
    NodeContainerConstIterator pointsEnd = this->m_TrialNodes.end();

    while( pointsIter != pointsEnd )
      {
      // get node from trial points container
      idx = pointsIter->first;

      // check if node index is within the output level set
      //if ( m_BufferedRegion.IsInside( idx ) )
        {
        // make this an initial trial point
        m_LabelImage->SetPixel( idx, Superclass::InitialTrialPoint );

        output->SetPixel(idx, pointsIter->second);

        this->m_Heap->Push( PriorityQueueElementType( idx, pointsIter->second ) );
        //m_TrialHeap.push(node);
        }
      ++pointsIter;
      }
    }
  }
// -----------------------------------------------------------------------------
}
#endif // __itkFastMarchingImageFilterBase_txx
