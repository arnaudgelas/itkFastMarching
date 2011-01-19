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

#ifndef __itkFastMarchingImageFilterBase_h
#define __itkFastMarchingImageFilterBase_h

#include "itkFastMarchingBase.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{

template< unsigned int VDimension, typename TInputPixel, typename TOutputPixel >
class FastMarchingImageFilterBase :
    public FastMarchingBase<
      ImageFastMarchingTraits< VDimension, TInputPixel, TOutputPixel >
    >
  {
public:
  typedef ImageFastMarchingTraits< VDimension, TInputPixel, TOutputPixel > Traits;

  typedef FastMarchingImageFilterBase            Self;
  typedef FastMarchingBase< Traits >             Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;


  typedef typename Superclass::InputDomainType     InputImageType;
  typedef typename Superclass::InputDomainPointer  InputImagePointer;
  typedef typename Superclass::InputPixelType      InputPixelType;

  typedef typename Superclass::OutputDomainType     OutputImageType;
  typedef typename Superclass::OutputDomainPointer  OutputImagePointer;
  typedef typename Superclass::OutputPixelType      OutputPixelType;
  typedef typename OutputImageType::SpacingType     OutputSpacingType;
  typedef typename OutputImageType::SizeType        OutputSizeType;

  typedef typename Superclass::NodeType NodeType;
  typedef typename Superclass::ElementIdentifier ElementIdentifier;

  typedef typename Superclass::PriorityQueueElementType PriorityQueueElementType;

  typedef typename Superclass::PriorityQueueType PriorityQueueType;
  typedef typename Superclass::PriorityQueuePointer PriorityQueuePointer;

  typedef typename Superclass::LabelType LabelType;

  itkStaticConstMacro( ImageDimension, unsigned int, VDimension );

  typedef Image< unsigned char, ImageDimension >  LabelImageType;
  typedef typename LabelImageType::Pointer        LabelImagePointer;

  typedef typename Superclass::NodeContainerType NodeContainerType;
  typedef typename Superclass::NodeContainerConstIterator NodeContainerConstIterator;

  /** \brief Set a binary mask to limit the propagation of front.
    *
    * Only location with null pixel values are used. Corresponding labels are set
    * to Forbidden to limit the propagation of the front and output values are set
    * to 0.
    *
    * \param[in] iImage input image
    * \tparam TPixel PixelType for the input image
   */
  template< typename TPixel >
  void SetBinaryMask( Image< TPixel, ImageDimension >* iImage )
    {
    OutputSizeType nullsize;
    nullsize.Fill( 0 );

    /*
    if( m_BufferedRegion.GetSize() == nullsize )
      {
      itkGenericExceptionMacro( << "m_BufferedRegion has not been set yet" );
      }
    else */
      {
      typedef Image< TPixel, ImageDimension > InternalImageType;
      typedef ImageRegionConstIteratorWithIndex< InternalImageType >
        InternalRegionIterator;
      InternalRegionIterator b_it( iImage, iImage->GetLargestPossibleRegion() );
      b_it.GoToBegin();

      OutputImageType* output = this->GetOutput();

      TPixel zero_value = NumericTraits< TPixel >::Zero;
      NodeType idx;

      while( !b_it.IsAtEnd() )
        {
        if( b_it.Get() == zero_value )
          {
          idx = b_it.GetIndex();
          //if ( m_BufferedRegion.IsInside( idx ) )
            {
            m_LabelImage->SetPixel(idx, Superclass::Forbidden );
            output->SetPixel (idx, zero_value );
            }
          }
        ++b_it;
        }
      this->Modified();
      }
    }

protected:
  FastMarchingImageFilterBase();
  ~FastMarchingImageFilterBase();

  struct InternalNodeStructure;

  NodeType m_StartIndex;
  NodeType m_LastIndex;

  LabelImagePointer m_LabelImage;

  LabelType GetLabelValueForGivenNode( NodeType iNode );
  void SetLabelValueForGivenNode( NodeType iNode, LabelType iLabel );
  void UpdateNeighbors( NodeType iNode );
  void UpdateValue( NodeType iValue );
  void CheckTopology( NodeType iNode );
  void InitializeOutput();
  double Solve( std::vector< InternalNodeStructure > iNeighbors );

private:
  FastMarchingImageFilterBase( const Self& );
  void operator = ( const Self& );
  };
}

#include "itkFastMarchingImageFilterBase.txx"
#endif // __itkFastMarchingImageFilterBase_h
