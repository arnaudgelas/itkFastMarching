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
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

namespace itk
{
// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
class
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
    InternalNodeStructure
  {
public:
    InternalNodeStructure( ) :
      m_Value( NumericTraits< TOutputPixel >::max() ) {}

    NodeType        m_Node;
    OutputPixelType m_Value;
    unsigned int    m_Axis;

    bool operator< ( const InternalNodeStructure& iRight ) const
      {
      return m_Value < iRight.m_Value;
      }
  };


// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
FastMarchingImageFilterBase()
  {
  OutputSizeType outputSize;
  outputSize.Fill(16);

  NodeType outputIndex;
  outputIndex.Fill(0);

  m_OutputRegion.SetSize(outputSize);
  m_OutputRegion.SetIndex(outputIndex);

  m_OutputOrigin.Fill(0.0);
  m_OutputSpacing.Fill(1.0);
  m_OutputDirection.SetIdentity();
  m_OverrideOutputInformation = false;

  m_LabelImage = LabelImageType::New();
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
char
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

  double solution = Solve( iNode, NodesUsed );

  if ( solution < this->m_LargeValue )
    {
    // write solution to m_OutputLevelSet
    OutputPixelType outputPixel = static_cast< OutputPixelType >( solution );
    this->GetOutput()->SetPixel(iNode, outputPixel);

    // insert point into trial heap
    m_LabelImage->SetPixel( iNode, Superclass::Trial );
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
Solve( NodeType iNode, std::vector< InternalNodeStructure > iNeighbors )
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
      static_cast< double >( this->GetInput()->GetPixel(iNode) ) /
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
bool
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
CheckTopology( NodeType iNode )
  {
  if( this->m_TopologyCheck != Superclass::None )
    {
    if( ( ImageDimension == 2 ) || ( ImageDimension == 3 ) )
      {
      bool wellComposednessViolation
        = this->DoesVoxelChangeViolateWellComposedness( iNode );

      bool strictTopologyViolation
         = this->DoesVoxelChangeViolateStrictTopology( iNode );

      if( ( this->m_TopologyCheck == Superclass::Strict ) &&
          ( wellComposednessViolation || strictTopologyViolation ) )
        {
        this->GetOutput()->SetPixel( iNode, this->m_TopologyValue );
        this->m_LabelImage->SetPixel( iNode, Superclass::Topology );
        return false;
        }
      if( this->m_TopologyCheck == Superclass::NoHandles )
        {
        if( wellComposednessViolation )
          {
          this->GetOutput()->SetPixel( iNode, this->m_TopologyValue );
          m_LabelImage->SetPixel( iNode, Superclass::Topology );
          return false;
          }
        if( strictTopologyViolation )
          {
          // check for handles
          typename NeighborhoodIteratorType::RadiusType radius;
          radius.Fill( 1 );
          NeighborhoodIteratorType ItL( radius, this->m_LabelImage,
            this->m_LabelImage->GetBufferedRegion() );
          ItL.SetLocation( iNode );

          NeighborhoodIterator<ConnectedComponentImageType> ItC(
                radius, this->m_ConnectedComponentImage,
                this->m_ConnectedComponentImage->GetBufferedRegion() );
          ItC.SetLocation( iNode );

          typename ConnectedComponentImageType::PixelType minLabel
              = NumericTraits<typename ConnectedComponentImageType::PixelType>::Zero;
          typename ConnectedComponentImageType::PixelType otherLabel
              = NumericTraits<typename ConnectedComponentImageType::PixelType>::Zero;

          bool doesChangeCreateHandle = false;

          for( unsigned int d = 0; d < ImageDimension; d++ )
            {
            if( ItL.GetNext( d ) == Superclass::Alive &&
                ItL.GetPrevious( d ) == Superclass::Alive )
              {
              if( ItC.GetNext( d ) == ItC.GetPrevious( d ) )
                {
                doesChangeCreateHandle = true;
                }
              else
                {
                minLabel = vnl_math_min( ItC.GetNext( d ), ItC.GetPrevious( d ) );
                otherLabel = vnl_math_max( ItC.GetNext( d ), ItC.GetPrevious( d ) );
                }
              break;
              }
            }
          if( doesChangeCreateHandle )
            {
            this->GetOutput()->SetPixel( iNode, this->m_TopologyValue );
            this->m_LabelImage->SetPixel( iNode, Superclass::Topology );
            return false;
            }
          else
            {
            for( ItC.GoToBegin(); !ItC.IsAtEnd(); ++ItC )
              {
              if( ItC.GetCenterPixel() == otherLabel )
                {
                ItC.SetCenterPixel( minLabel );
                }
              }
            }
          }
        }
      }
    else
      {
      itkWarningMacro( << "CheckTopology has not be implemented for Dimension != 2 and != 3."
                    << "m_TopologyCheck should be set to None." );
      }
    }
  return true;
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
  m_BufferedRegion = output->GetBufferedRegion();
  m_StartIndex = m_BufferedRegion.GetIndex();
  m_LastIndex = m_StartIndex + m_BufferedRegion.GetSize();

  typename OutputImageType::OffsetType offset;
  offset.Fill(1);
  m_LastIndex -= offset;

  // Checking for handles only requires an image to keep track of
  // connected components.
  if( this->m_TopologyCheck == Superclass::NoHandles )
    {
    m_ConnectedComponentImage = ConnectedComponentImageType::New();
    m_ConnectedComponentImage->SetOrigin( output->GetOrigin() );
    m_ConnectedComponentImage->SetSpacing( output->GetSpacing() );
    m_ConnectedComponentImage->SetRegions( output->GetBufferedRegion() );
    m_ConnectedComponentImage->SetDirection( output->GetDirection() );
    m_ConnectedComponentImage->Allocate();
    m_ConnectedComponentImage->FillBuffer( 0 );
    }

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
      if ( m_BufferedRegion.IsInside( idx ) )
        {
        // make this an alive point
        m_LabelImage->SetPixel(idx, Superclass::Alive );

        if( this->m_TopologyCheck == Superclass::NoHandles )
          {
          m_ConnectedComponentImage->SetPixel( idx,
            NumericTraits<unsigned int>::One );
          }

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
      if ( m_BufferedRegion.IsInside( idx ) )
        {
        // make this an alive point
        m_LabelImage->SetPixel(idx, Superclass::Forbidden );
        output->SetPixel (idx, zero );
        }

      ++p_it;
      }
    }

  if( this->m_TopologyCheck == Superclass::NoHandles )
    {
    // Now create the connected component image and relabel such that labels
    // are 1, 2, 3, ...
    typedef ConnectedComponentImageFilter<ConnectedComponentImageType,
      ConnectedComponentImageType> ConnectedComponentFilterType;
    typename ConnectedComponentFilterType::Pointer connecter
        = ConnectedComponentFilterType::New();
    connecter->SetInput( m_ConnectedComponentImage );

    typedef RelabelComponentImageFilter<ConnectedComponentImageType,
        ConnectedComponentImageType> RelabelerType;
    typename RelabelerType::Pointer relabeler = RelabelerType::New();
    relabeler->SetInput( connecter->GetOutput() );
    relabeler->Update();

    this->m_ConnectedComponentImage = relabeler->GetOutput();
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
      if ( m_BufferedRegion.IsInside( idx ) )
        {
        // make this an initial trial point
        m_LabelImage->SetPixel( idx, Superclass::InitialTrial );

        output->SetPixel(idx, pointsIter->second);

        //this->m_Heap->Push( PriorityQueueElementType( idx, pointsIter->second ) );
        this->m_Heap.push( NodePairType( idx, pointsIter->second ) );
        }
      ++pointsIter;
      }
    }
  // initialize indices if this->m_TopologyCheck is activated
  if( this->m_TopologyCheck != Superclass::None )
    {
    if( ImageDimension == 2 )
      {
      InitializeIndices2D();
      }
    else
      {
      if( ImageDimension == 3 )
        {
        InitializeIndices3D();
        }
      else
        {
        itkWarningMacro(
              << "Topology checking is only valid for level set dimensions of 2 and 3" );
        }
      }
    }
  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
bool
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
DoesVoxelChangeViolateWellComposedness( NodeType idx )
{
  bool isChangeWellComposed = false;
  if( ImageDimension == 2 )
    {
    isChangeWellComposed = this->IsChangeWellComposed2D( idx );
    }
  else  // ImageDimension == 3
    {
    isChangeWellComposed = this->IsChangeWellComposed3D( idx );
    }

  return !isChangeWellComposed;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
bool
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
DoesVoxelChangeViolateStrictTopology( NodeType idx )
{
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  NeighborhoodIteratorType It( radius, this->m_LabelImage,
    this->m_LabelImage->GetBufferedRegion() );
  It.SetLocation( idx );

  unsigned int numberOfCriticalC3Configurations = 0;
  unsigned int numberOfFaces = 0;

  for( unsigned int d = 0; d < ImageDimension; d++ )
    {
    if( It.GetNext( d ) == Superclass::Alive )
      {
      ++numberOfFaces;
      }
    if( It.GetPrevious( d ) == Superclass::Alive )
      {
      ++numberOfFaces;
      }
    if( It.GetNext( d ) == Superclass::Alive &&
        It.GetPrevious( d ) == Superclass::Alive )
      {
      ++numberOfCriticalC3Configurations;
      }
    }

  if( ( numberOfCriticalC3Configurations > 0 ) &&
      ( numberOfFaces % 2 == 0 ) &&
      ( numberOfCriticalC3Configurations * 2 == numberOfFaces ) )
    {
    return true;
    }
  return false;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
bool
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
IsChangeWellComposed2D( NodeType idx )
{
  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );

  NeighborhoodIteratorType It( radius, this->m_LabelImage,
    this->m_LabelImage->GetBufferedRegion() );
  It.SetLocation( idx );

  std::vector<bool> neighborhoodPixels( 9 );

  // Check for critical configurations: 4 90-degree rotations

  for ( unsigned int i = 0; i < 4; i++ )
    {
    for ( unsigned int j = 0; j < 9; j++ )
      {
      neighborhoodPixels[j] =
        ( It.GetPixel( this->m_RotationIndices[i][j] ) != Superclass::Alive );
      if( this->m_RotationIndices[i][j] == 4 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }

    if( this->IsCriticalC1Configuration2D( neighborhoodPixels )
      || this->IsCriticalC2Configuration2D( neighborhoodPixels )
      || this->IsCriticalC3Configuration2D( neighborhoodPixels )
      || this->IsCriticalC4Configuration2D( neighborhoodPixels ) )
      {
      return false;
      }
    }

  // Check for critical configurations: 2 reflections
  //  Note that the reflections for the C1 and C2 cases
  //  are covered by the rotation cases above (except
  //  in the case of FullInvariance == false.

  for ( unsigned int i = 0; i < 2; i++ )
    {
    for ( unsigned int j = 0; j < 9; j++ )
      {
      neighborhoodPixels[j] =
        ( It.GetPixel( this->m_ReflectionIndices[i][j] ) != Superclass::Alive );
      if( this->m_ReflectionIndices[i][j] == 4 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
//    if( !this->m_FullInvariance
//      && ( this->IsCriticalC1Configuration2D( neighborhoodPixels )
//        || this->IsCriticalC2Configuration2D( neighborhoodPixels ) ) )
//      {
//      return false;
//      }
    if( this->IsCriticalC3Configuration2D( neighborhoodPixels )
      || this->IsCriticalC4Configuration2D( neighborhoodPixels ) )
      {
      return false;
      }
    }
  return true;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
bool
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
IsCriticalC1Configuration2D( std::vector<bool> neighborhood )
{
  return ( !neighborhood[0] &&  neighborhood[1] &&
            neighborhood[3] && !neighborhood[4] &&
           !neighborhood[8] );
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
bool
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
IsCriticalC2Configuration2D( std::vector<bool> neighborhood )
{
  return ( !neighborhood[0] &&  neighborhood[1] &&
            neighborhood[3] && !neighborhood[4] &&
            neighborhood[8] &&
           ( neighborhood[5] || neighborhood[7] ) );
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
bool
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
IsCriticalC3Configuration2D( std::vector<bool> neighborhood )
{
  return ( !neighborhood[0] &&  neighborhood[1] &&
            neighborhood[3] && !neighborhood[4] &&
           !neighborhood[5] &&  neighborhood[6] &&
           !neighborhood[7] &&  neighborhood[8] );
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
bool
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
IsCriticalC4Configuration2D( std::vector<bool> neighborhood )
{
  return ( !neighborhood[0] &&  neighborhood[1] &&
            neighborhood[3] && !neighborhood[4] &&
           !neighborhood[5] && !neighborhood[6] &&
           !neighborhood[7] &&  neighborhood[8] );
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
void
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
InitializeIndices2D()
{
  this->m_RotationIndices[0].SetSize( 9 );
  this->m_RotationIndices[1].SetSize( 9 );
  this->m_RotationIndices[2].SetSize( 9 );
  this->m_RotationIndices[3].SetSize( 9 );

  this->m_RotationIndices[0][0] = 0;
  this->m_RotationIndices[0][1] = 1;
  this->m_RotationIndices[0][2] = 2;
  this->m_RotationIndices[0][3] = 3;
  this->m_RotationIndices[0][4] = 4;
  this->m_RotationIndices[0][5] = 5;
  this->m_RotationIndices[0][6] = 6;
  this->m_RotationIndices[0][7] = 7;
  this->m_RotationIndices[0][8] = 8;

  this->m_RotationIndices[1][0] = 2;
  this->m_RotationIndices[1][1] = 5;
  this->m_RotationIndices[1][2] = 8;
  this->m_RotationIndices[1][3] = 1;
  this->m_RotationIndices[1][4] = 4;
  this->m_RotationIndices[1][5] = 7;
  this->m_RotationIndices[1][6] = 0;
  this->m_RotationIndices[1][7] = 3;
  this->m_RotationIndices[1][8] = 6;

  this->m_RotationIndices[2][0] = 8;
  this->m_RotationIndices[2][1] = 7;
  this->m_RotationIndices[2][2] = 6;
  this->m_RotationIndices[2][3] = 5;
  this->m_RotationIndices[2][4] = 4;
  this->m_RotationIndices[2][5] = 3;
  this->m_RotationIndices[2][6] = 2;
  this->m_RotationIndices[2][7] = 1;
  this->m_RotationIndices[2][8] = 0;

  this->m_RotationIndices[3][0] = 6;
  this->m_RotationIndices[3][1] = 3;
  this->m_RotationIndices[3][2] = 0;
  this->m_RotationIndices[3][3] = 7;
  this->m_RotationIndices[3][4] = 4;
  this->m_RotationIndices[3][5] = 1;
  this->m_RotationIndices[3][6] = 8;
  this->m_RotationIndices[3][7] = 5;
  this->m_RotationIndices[3][8] = 2;

  this->m_ReflectionIndices[0].SetSize( 9 );
  this->m_ReflectionIndices[1].SetSize( 9 );

  this->m_ReflectionIndices[0][0] = 6;
  this->m_ReflectionIndices[0][1] = 7;
  this->m_ReflectionIndices[0][2] = 8;
  this->m_ReflectionIndices[0][3] = 3;
  this->m_ReflectionIndices[0][4] = 4;
  this->m_ReflectionIndices[0][5] = 5;
  this->m_ReflectionIndices[0][6] = 0;
  this->m_ReflectionIndices[0][7] = 1;
  this->m_ReflectionIndices[0][8] = 2;

  this->m_ReflectionIndices[1][0] = 2;
  this->m_ReflectionIndices[1][1] = 1;
  this->m_ReflectionIndices[1][2] = 0;
  this->m_ReflectionIndices[1][3] = 5;
  this->m_ReflectionIndices[1][4] = 4;
  this->m_ReflectionIndices[1][5] = 3;
  this->m_ReflectionIndices[1][6] = 8;
  this->m_ReflectionIndices[1][7] = 7;
  this->m_ReflectionIndices[1][8] = 6;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
bool
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
IsChangeWellComposed3D( NodeType idx )
{
  std::vector<bool> neighborhoodPixels( 8 );

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill( 1 );
  NeighborhoodIteratorType It( radius, this->m_LabelImage,
    this->m_LabelImage->GetRequestedRegion() );
  It.SetLocation( idx );

  // Check for C1 critical configurations
  for ( unsigned int i = 0; i < 12; i++ )
    {
    for ( unsigned int j = 0; j < 4; j++ )
      {
      neighborhoodPixels[j]
        = ( It.GetPixel( m_C1Indices[i][j] ) == Superclass::Alive );
      if( m_C1Indices[i][j] == 13 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
    if( this->IsCriticalC1Configuration3D( neighborhoodPixels ) )
      {
      return false;
      }
    }

  // Check for C2 critical configurations
  for ( unsigned int i = 0; i < 8; i++ )
    {
    for ( unsigned int j = 0; j < 8; j++ )
      {
      neighborhoodPixels[j]
        = ( It.GetPixel( m_C2Indices[i][j] ) == Superclass::Alive );
      if( m_C2Indices[i][j] == 13 )
        {
        neighborhoodPixels[j] = !neighborhoodPixels[j];
        }
      }
    if( IsCriticalC2Configuration3D( neighborhoodPixels ) )
      {
      return false;
      }
    }

  return true;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
bool
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
IsCriticalC1Configuration3D( std::vector<bool> neighborhood )
{
  return ( (  neighborhood[0] &&  neighborhood[1] &&
             !neighborhood[2] && !neighborhood[3] ) ||
           ( !neighborhood[0] && !neighborhood[1] &&
              neighborhood[2] &&  neighborhood[3] ) );
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
unsigned int
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
IsCriticalC2Configuration3D( std::vector<bool> neighborhood )
{
  // Check if Type 1 or Type 2
  for ( unsigned int i = 0; i < 4; i++ )
    {
    bool isC2 = false;
    if( neighborhood[2*i] == neighborhood[2*i+1] )
      {
      isC2 = true;
      for ( unsigned int j = 0; j < 8; j++ )
        {
        if( neighborhood[j] == neighborhood[2*i] &&
               j != 2*i && j != 2*i+1 )
          {
          isC2 = false;
          }
        }
      }
    if( isC2 )
      {
      if( neighborhood[2*i] )
        {
        return 1;
        }
      else
        {
        return 2;
        }
      }
    }

  return 0;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
void
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
InitializeIndices3D()
{
  for ( unsigned int i = 0; i <  12; i++ )
    {
    this->m_C1Indices[i].SetSize( 4 );
    }
  for ( unsigned int i = 0; i <  8; i++ )
    {
    this->m_C2Indices[i].SetSize( 8 );
    }

  this->m_C1Indices[0][0] = 1;
  this->m_C1Indices[0][1] = 13;
  this->m_C1Indices[0][2] = 4;
  this->m_C1Indices[0][3] = 10;

  this->m_C1Indices[1][0] = 9;
  this->m_C1Indices[1][1] = 13;
  this->m_C1Indices[1][2] = 10;
  this->m_C1Indices[1][3] = 12;

  this->m_C1Indices[2][0] = 3;
  this->m_C1Indices[2][1] = 13;
  this->m_C1Indices[2][2] = 4;
  this->m_C1Indices[2][3] = 12;

  this->m_C1Indices[3][0] = 4;
  this->m_C1Indices[3][1] = 14;
  this->m_C1Indices[3][2] = 5;
  this->m_C1Indices[3][3] = 13;

  this->m_C1Indices[4][0] = 12;
  this->m_C1Indices[4][1] = 22;
  this->m_C1Indices[4][2] = 13;
  this->m_C1Indices[4][3] = 21;

  this->m_C1Indices[5][0] = 13;
  this->m_C1Indices[5][1] = 23;
  this->m_C1Indices[5][2] = 14;
  this->m_C1Indices[5][3] = 22;

  this->m_C1Indices[6][0] = 4;
  this->m_C1Indices[6][1] = 16;
  this->m_C1Indices[6][2] = 7;
  this->m_C1Indices[6][3] = 13;

  this->m_C1Indices[7][0] = 13;
  this->m_C1Indices[7][1] = 25;
  this->m_C1Indices[7][2] = 16;
  this->m_C1Indices[7][3] = 22;

  this->m_C1Indices[8][0] = 10;
  this->m_C1Indices[8][1] = 22;
  this->m_C1Indices[8][2] = 13;
  this->m_C1Indices[8][3] = 19;

  this->m_C1Indices[9][0] = 12;
  this->m_C1Indices[9][1] = 16;
  this->m_C1Indices[9][2] = 13;
  this->m_C1Indices[9][3] = 15;

  this->m_C1Indices[10][0] = 13;
  this->m_C1Indices[10][1] = 17;
  this->m_C1Indices[10][2] = 14;
  this->m_C1Indices[10][3] = 16;

  this->m_C1Indices[11][0] = 10;
  this->m_C1Indices[11][1] = 14;
  this->m_C1Indices[11][2] = 11;
  this->m_C1Indices[11][3] = 13;

  this->m_C2Indices[0][0] = 0;
  this->m_C2Indices[0][1] = 13;
  this->m_C2Indices[0][2] = 1;
  this->m_C2Indices[0][3] = 12;
  this->m_C2Indices[0][4] = 3;
  this->m_C2Indices[0][5] = 10;
  this->m_C2Indices[0][6] = 4;
  this->m_C2Indices[0][7] = 9;

  this->m_C2Indices[4][0] = 9;
  this->m_C2Indices[4][1] = 22;
  this->m_C2Indices[4][2] = 10;
  this->m_C2Indices[4][3] = 21;
  this->m_C2Indices[4][4] = 12;
  this->m_C2Indices[4][5] = 19;
  this->m_C2Indices[4][6] = 13;
  this->m_C2Indices[4][7] = 18;

  for ( unsigned int i = 1; i < 4; i++ )
    {
    int addend;
    if( i == 2 )
      {
      addend = 2;
      }
    else
      {
      addend = 1;
      }
    for ( unsigned int j = 0; j < 8; j++ )
      {
      this->m_C2Indices[i  ][j] = this->m_C2Indices[i-1][j] + addend;
      this->m_C2Indices[i+4][j] = this->m_C2Indices[i+3][j] + addend;
      }
    }
}
// -----------------------------------------------------------------------------
}
#endif // __itkFastMarchingImageFilterBase_txx
