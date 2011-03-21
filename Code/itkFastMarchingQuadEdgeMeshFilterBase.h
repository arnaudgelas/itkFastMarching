/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009-11

 Copyright (c) 2009-11, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#ifndef __itkFastMarchingQuadEdgeMeshFilterBase_h
#define __itkFastMarchingQuadEdgeMeshFilterBase_h

#include "itkFastMarchingBase.h"
#include "itkFastMarchingTraits.h"

namespace itk
{

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
class FastMarchingQuadEdgeMeshFilterBase :
    public FastMarchingBase<
      MeshFastMarchingTraits< VDimension,
        TInputPixel,
        TOutputPixel,
        TOutputMeshTraits > >
{
public:

  typedef MeshFastMarchingTraits< VDimension,
    TInputPixel,
    TOutputPixel,
    TOutputMeshTraits > Traits;

  typedef FastMarchingQuadEdgeMeshFilterBase     Self;
  typedef FastMarchingBase< Traits >             Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FastMarchingQuadEdgeMeshFilterBase, FastMarchingBase);

  typedef typename Superclass::InputDomainType     InputImageType;
  typedef typename Superclass::InputDomainPointer  InputImagePointer;
  typedef typename Superclass::InputPixelType      InputPixelType;

  typedef typename Superclass::OutputDomainType     OutputImageType;
  typedef typename Superclass::OutputDomainPointer  OutputImagePointer;
  typedef typename Superclass::OutputPixelType      OutputPixelType;

  typedef std::map< NodeType, LabelType > NodeLabelMapType;
  typedef typename NodeLabelMapType::iterator NodeLabelMapIterator;
  typedef typename NodeLabelMapType::const_iterator NodeLabelMapConstIterator;

protected:
  FastMarchingQuadEdgeMeshFilterBase() {}
  virtual ~FastMarchingQuadEdgeMeshFilterBase() {}


  IdentifierType GetTotalNumberOfNodes() const
    {
    this->GetInput()->GetNumberOfPoints();
    }

  OutputPixelType GetOutputValue( OutputDomainType* oDomain,
                                  const NodeType& iNode ) const
    {
    return oDomain->GetPointData( iNode );
    }

  char GetLabelValueForGivenNode( const NodeType& iNode ) const
    {
    NodeLabelMapIterator it = m_Label.find( iNode );

    if( it != m_Label.end() )
      {
      return iLabel;
      }
    else
      {
      return Superclass::Far;
      }
    }

  void SetLabelValueForGivenNode( const NodeType& iNode,
                                  const LabelType& iLabel )
    {
    m_Label[iNode] = iLabel;
    }

  void UpdateNeighbors( OutputDomainType* oDomain,
                        const NodeType& iNode )
    {

    }

  void UpdateValue( OutputDomainType* oDomain,
                    const NodeType& iNode )
    {

    }

  bool CheckTopology( OutputDomainType* oDomain,
                      const NodeType& iNode )
    {
    itkWarningMacro( << "Constrained topology on Mesh is not implemented" );
    }

  void InitializeOutput( OutputDomainType* oDomain )
    {
    this->CopyInputMeshToOutputMesh();

    OutputMeshPointer output = this->GetOutput();

      {
      OutputPointsContainerPointer points = output->GetPoints();
      OutputPointsContainerIterator p_it = points->Begin();
      OutputPointsContainerIterator p_end = points->End();

      while( p_it != p_end )
        {
        output->SetPointData( p_it->Index(), this->m_LargeValue );
        ++p_it;
        }
      }

    m_Label.clear();

    if ( this->m_AlivePoints )
      {
      NodePairContainerConstIterator pointsIter = this->m_AlivePoints->Begin();
      NodePairContainerConstIterator pointsEnd = this->m_AlivePoints->End();

      while( pointsIter != pointsEnd )
        {
        // get node from alive points container
        idx = pointsIter->Value().GetNode();
        outputPixel = pointsIter->Value().GetValue();

        m_Label[idx] = Superclass::Alive;
        output->SetPointData( idx, outputPixel );

        ++pointsIter;
        }
      }
    if( this->m_ForbiddenPoints )
      {
      NodeContainerConstIterator node_it = this->m_ForbiddenPoints->Begin();
      NodeContainerConstIterator node_end = this->m_ForbiddenPoints->End();

      OutputPixelType zero = NumericTraits< OutputPixelType >::Zero;

      while( node_it != node_end )
        {
        idx = node_it->Value();

        m_Label[idx] = Superclass::Forbidden;
        output->SetPointData( idx, zero );

        ++node_it;
        }
      }
    if ( this->m_TrialPoints )
      {
      NodePairContainerConstIterator pointsIter = this->m_TrialPoints->Begin();
      NodePairContainerConstIterator pointsEnd = this->m_TrialPoints->End();

      while( pointsIter != pointsEnd )
        {
        idx = pointsIter->Value().GetNode();
        outputPixel = pointsIter->Value().GetValue();

        m_Label[idx] = Superclass::InitialTrial;

        output->SetPointData( idx, outputPixel );

        this->m_Heap.push( pointsIter->Value() );

        ++pointsIter;
        }
      }

  };


private:
  FastMarchingQuadEdgeMeshFilterBase( const Self& );
  void operator = ( const Self& );
};
}
#endif // __itkFastMarchingQuadEdgeMeshFilterBase_h
