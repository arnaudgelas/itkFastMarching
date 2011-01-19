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

namespace itk
{
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
    // update left neighbor
    if ( iNode[j] > m_StartIndex[j] )
      {
      neighIndex[j] = iNode[j] - 1;
      }

    label = m_LabelImage->GetPixel(neighIndex);

    if ( ( label != Superclass::Alive ) &&
         ( label != Superclass::InitialTrial ) &&
         ( label != Superclass::Forbidden ) )
      {
      this->UpdateValue( neighIndex );
      }

    // update right neighbor
    if ( iNode[j] < m_LastIndex[j] )
      {
      neighIndex[j] = iNode[j] + 1;
      }

    label = m_LabelImage->GetPixel(neighIndex);

    if ( ( label != Superclass::Alive ) &&
         ( label != Superclass::InitialTrial ) &&
         ( label != Superclass::Outside ) )
      {
      this->UpdateValue( neighIndex );
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

  }
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
void
FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
CheckTopology( NodeType iNode )
  {}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
void FastMarchingImageFilterBase< VDimension, TInputPixel, TOutputPixel >::
InitializeOutput()
  {}
// -----------------------------------------------------------------------------
}
#endif // __itkFastMarchingImageFilterBase_txx
