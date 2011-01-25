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
#ifndef __itkFastMarchingExtensionImageFilterBase_txx
#define __itkFastMarchingExtensionImageFilterBase_txx

#include "itkFastMarchingExtensionImageFilterBase.h"

namespace itk
{
template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TCriterion,
         class TAuxValue,
         unsigned int VAuxDimension >
FastMarchingExtensionImageFilterBase< VDimension, TInputPixel, TOutputPixel,
  TCriterion, TAuxValue, VAuxDimension >
::FastMarchingExtensionImageFilterBase()
{
  this->ProcessObject::SetNumberOfRequiredOutputs(1 + AuxDimension);

  AuxImagePointer ptr;
  for ( unsigned int k = 0; k < VAuxDimension; k++ )
    {
    ptr = AuxImageType::New();
    this->ProcessObject::SetNthOutput( k + 1, ptr.GetPointer() );
    }
}

template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TCriterion,
         class TAuxValue,
         unsigned int VAuxDimension >
void
FastMarchingExtensionImageFilterBase< VDimension, TInputPixel, TOutputPixel,
  TCriterion, TAuxValue, VAuxDimension >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Aux alive values: ";
  //os << m_AuxAliveValues.GetPointer() << std::endl;
  os << indent << "Aux trail values: ";
  //os << m_AuxTrialValues.GetPointer() << std::endl;
}


 template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TCriterion,
         class TAuxValue,
         unsigned int VAuxDimension >
typename FastMarchingExtensionImageFilterBase< VDimension, TInputPixel, TOutputPixel,
  TCriterion, TAuxValue, VAuxDimension >
::AuxImageType *
FastMarchingExtensionImageFilterBase< VDimension, TInputPixel, TOutputPixel,
  TCriterion, TAuxValue, VAuxDimension >
::GetAuxiliaryImage(unsigned int idx)
{
  if ( idx >= AuxDimension || this->GetNumberOfOutputs() < idx + 2 )
    {
    return NULL;
    }

  return static_cast< AuxImageType * >( this->ProcessObject::GetOutput(idx + 1) );
}

/*
 *
 */
template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TCriterion,
         class TAuxValue,
         unsigned int VAuxDimension >
void
FastMarchingExtensionImageFilterBase< VDimension, TInputPixel, TOutputPixel,
  TCriterion, TAuxValue, VAuxDimension >
::GenerateOutputInformation()
{
  // call the superclass implementation of this function
  this->Superclass::GenerateOutputInformation();

  // set the size of all the auxiliary outputs
  // to be the same as the primary output
  OutputImageType* primaryOutput = this->GetOutput();
  for ( unsigned int k = 0; k < VAuxDimension; k++ )
    {
    AuxImageType *ptr = this->GetAuxiliaryImage(k);
    ptr->CopyInformation(primaryOutput);
    }
}

/*
 *
 */
template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TCriterion,
         class TAuxValue,
         unsigned int VAuxDimension >
void
FastMarchingExtensionImageFilterBase< VDimension, TInputPixel, TOutputPixel,
  TCriterion, TAuxValue, VAuxDimension >
::EnlargeOutputRequestedRegion(
  DataObject *itkNotUsed(output) )
{
  // This filter requires all of the output images in the buffer.
  for ( unsigned int j = 0; j < this->GetNumberOfOutputs(); j++ )
    {
    if ( this->ProcessObject::GetOutput(j) )
      {
      this->ProcessObject::GetOutput(j)->SetRequestedRegionToLargestPossibleRegion();
      }
    }
}

/*
 *
 */
template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TCriterion,
         class TAuxValue,
         unsigned int VAuxDimension >
void
FastMarchingExtensionImageFilterBase< VDimension, TInputPixel, TOutputPixel,
  TCriterion, TAuxValue, VAuxDimension >
::InitializeOutput(OutputImageType* oImage)
{
  this->Superclass::InitializeOutput( oImage );

  if ( m_AuxAliveValues.size() != ( this->m_AliveNodes.size() ) )
    {
    itkExceptionMacro(<< "in Initialize(): AuxAliveValues is the wrong size");
    }

  if ( m_AuxTrialValues.size() != this->m_TrialNodes.size() )
    {
    itkExceptionMacro(<< "in Initialize(): AuxTrialValues is the wrong size");
    }

  AuxImagePointer auxImages[AuxDimension];

  // allocate memory for the auxiliary outputs
  for ( unsigned int k = 0; k < VAuxDimension; k++ )
    {
    AuxImageType *ptr = this->GetAuxiliaryImage(k);
    ptr->SetBufferedRegion( ptr->GetRequestedRegion() );
    ptr->Allocate();
    auxImages[k] = ptr;
    }

  // set all alive points to alive
  typename Superclass::NodeType node;

  AuxValueVectorType auxVec;

  if ( !m_AuxAliveValues.empty() )
    {
    AuxValueContainerConstIterator auxIter = m_AuxAliveValues.begin();
    NodeContainerConstIterator pointsIter =  this->m_AliveNodes.begin();
    NodeContainerConstIterator pointsEnd =  this->m_AliveNodes.end();

    while( pointsIter != pointsEnd )
      {
      node = pointsIter->first;
      auxVec = *auxIter;

      // check if node index is within the output level set
      if ( this->m_BufferedRegion.IsInside( node ) )
        {
        for ( unsigned int k = 0; k < VAuxDimension; k++ )
          {
          auxImages[k]->SetPixel( node, auxVec[k] );
          }
        }
      ++pointsIter;
      ++auxIter;
      } // end container while
    }   // if AuxAliveValues set

  if ( !m_AuxTrialValues.empty() )
    {
    AuxValueContainerConstIterator auxIter = m_AuxTrialValues.begin();
    NodeContainerConstIterator pointsIter =  this->m_TrialNodes.begin();
    NodeContainerConstIterator pointsEnd =  this->m_TrialNodes.end();

    while ( pointsIter != pointsEnd )
      {
      node = pointsIter->first;
      auxVec = *auxIter;

      // check if node index is within the output level set
      if ( this->m_BufferedRegion.IsInside( node ) )
        {
        for ( unsigned int k = 0; k < VAuxDimension; k++ )
          {
          auxImages[k]->SetPixel( node, auxVec[k]);
          }
        }
      ++pointsIter;
      ++auxIter;
      } // end container loop
    }   // if AuxTrialValues set
}

template< unsigned int VDimension,
         typename TInputPixel,
         typename TOutputPixel,
         class TCriterion,
         class TAuxValue,
         unsigned int VAuxDimension >
void
FastMarchingExtensionImageFilterBase< VDimension, TInputPixel, TOutputPixel,
  TCriterion, TAuxValue, VAuxDimension >
::UpdateValue( OutputImageType* oImage, const NodeType& iNode )
{
  // A extension value at node is choosen such that
  // grad(F) dot_product grad(Phi) = 0
  // where F is the extended speed function and Phi is
  // the level set function.
  //
  // The extension value can approximated as a weighted
  // sum of the values from nodes used in the calculation
  // of the distance by the superclass.
  //
  // For more detail see Chapter 11 of
  // "Level Set Methods and Fast Marching Methods", J.A. Sethian,
  // Cambridge Press, Second edition, 1999.

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
      if ( ( neighbor_node[j] > this->m_LastIndex[j] ) ||
           ( neighbor_node[j] < this->m_StartIndex[j] ) )
        {
        continue;
        }

      if ( this->m_LabelImage->GetPixel( neighbor_node ) == Superclass::Alive )
        {
        neighValue =
            static_cast< OutputPixelType >( oImage->GetPixel( neighbor_node ) );

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

  OutputPixelType outputPixel =
      static_cast< OutputPixelType >( Solve( oImage, iNode, NodesUsed ) );

  if ( outputPixel < this->m_LargeValue )
    {
    oImage->SetPixel(iNode, outputPixel);

    // insert point into trial heap
    this->m_LabelImage->SetPixel( iNode, Superclass::Trial );

    //node.SetValue( outputPixel );
    //node.SetIndex( index );
    //m_TrialHeap.push(node);
    this->m_Heap.push( NodePairType( iNode, outputPixel ) );

    // update auxiliary values
    for ( unsigned int k = 0; k < VAuxDimension; k++ )
      {
      double       numer = 0.0;
      double       denom = 0.;
      AuxValueType auxVal;

      for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
        temp_node = NodesUsed[j];

        if ( outputPixel < temp_node.m_Value )
          {
          break;
          }

        auxVal = this->GetAuxiliaryImage(k)->GetPixel( temp_node.m_Node );
        numer += auxVal * ( outputPixel - temp_node.m_Value );
        denom += outputPixel - temp_node.m_Value;
        }

      if ( denom > 0 )
        {
        auxVal = static_cast< AuxValueType >( numer / denom );
        }
      else
        {
        auxVal = NumericTraits< AuxValueType >::Zero;
        }

      this->GetAuxiliaryImage(k)->SetPixel(iNode, auxVal);
      }
    }
}
} // namespace itk

#endif
