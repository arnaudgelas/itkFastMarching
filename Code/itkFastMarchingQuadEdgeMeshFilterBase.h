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

#include "vnl/vnl_math.h"

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
        TInputMeshTraits,
        TOutputPixel,
        TOutputMeshTraits > >
{
public:

  typedef MeshFastMarchingTraits< VDimension,
    TInputPixel,
    TInputMeshTraits,
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

  typedef typename Superclass::InputDomainType     InputMeshType;
  typedef typename Superclass::InputDomainPointer  InputMeshPointer;
  typedef typename Superclass::InputPixelType      InputPixelType;
  typedef typename InputMeshType::PointType        InputPointType;
  typedef typename InputMeshType::PointIdentifier  InputPointIdentifierType;

  typedef typename Superclass::OutputDomainType     OutputMeshType;
  typedef typename Superclass::OutputDomainPointer  OutputMeshPointer;
  typedef typename Superclass::OutputPixelType      OutputPixelType;
  typedef typename OutputMeshType::PointType        OutputPointType;
  typedef typename OutputPointType::VectorType      OutputVectorType;
  typedef typename OutputVectorType::RealValueType  OutputVectorRealType;
  typedef typename OutputMeshType::QEType           OutputQEType;
  typedef typename OutputMeshType::PointIdentifier  OutputPointIdentifierType;
  typedef typename OutputMeshType::PointsContainer  OutputPointsContainer;
  typedef typename OutputPointsContainer::Pointer   OutputPointsContainerPointer;
  typedef typename OutputPointsContainer::Iterator  OutputPointsContainerIterator;


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

  typedef typename Superclass::LabelType LabelType;

  typedef std::map< NodeType, LabelType > NodeLabelMapType;
  typedef typename NodeLabelMapType::iterator NodeLabelMapIterator;
  typedef typename NodeLabelMapType::const_iterator NodeLabelMapConstIterator;

protected:

  FastMarchingQuadEdgeMeshFilterBase() {}
  virtual ~FastMarchingQuadEdgeMeshFilterBase() {}

  NodeLabelMapType m_Label;

  IdentifierType GetTotalNumberOfNodes() const
    {
    return this->GetInput()->GetNumberOfPoints();
    }

  OutputPixelType GetOutputValue( OutputMeshType* oDomain,
                                  const NodeType& iNode ) const
    {
    OutputPixelType outputValue;
    oDomain->GetPointData( iNode, &outputValue );
    return outputValue;
    }

  char GetLabelValueForGivenNode( const NodeType& iNode ) const
    {
    NodeLabelMapConstIterator it = m_Label.find( iNode );

    if( it != m_Label.end() )
      {
      return it->second;
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

  void UpdateNeighbors( OutputMeshType* oMesh,
                        const NodeType& iNode )
    {
    OutputPointType p;
    oMesh->GetPoint( iNode, &p );

    OutputQEType* qe = p.GetEdge();

    if( qe )
      {
      OutputQEType *qe_it = qe;
      qe_it = qe_it->GetOnext();

      while( qe_it != qe )
        {
        if( qe_it )
          {
          OutputPointIdentifierType neigh_id = qe_it->GetDestination();

          char label = this->GetLabelValueForGivenNode( neigh_id );

          if ( ( label != Superclass::Alive ) &&
               ( label != Superclass::InitialTrial ) &&
               ( label != Superclass::Forbidden ) )
            {
            this->UpdateValue( oMesh, neigh_id );
            }
          }
        qe_it = qe_it->GetOnext();
        }
      }
    else
      {
      // throw one exception here
      }
    }



  void UpdateValue( OutputMeshType* oMesh,
                    const NodeType& iNode )
    {
    OutputPointType p;
    oMesh->GetPoint( iNode, &p );

    OutputQEType* qe = p.GetEdge();

    OutputPixelType outputPixel = this->m_LargeValue;

    if( qe )
      {
      OutputQEType *qe_it = qe;
      qe_it = qe_it->GetOnext();

      do
        {
        OutputPointIdentifierType id1 = qe_it->GetDestination();

        OutputQEType *qe_it2 = qe_it->GetOnext();

        if( qe_it2 )
          {
          OutputPointIdentifierType id2 = qe_it2->GetDestination();

          char label1 = this->GetLabelValueForGivenNode( id1 );
          char label2 = this->GetLabelValueForGivenNode( id2 );

          if( ( label1 != Superclass::Far ) &&
              ( label2 != Superclass::Far ) )
            {
            OutputPointType q1 = oMesh->GetPoint( id1 );
            OutputPointType q2 = oMesh->GetPoint( id2 );

            OutputPixelType val1 = this->GetOutputValue( oMesh, id1 );
            OutputPixelType val2 = this->GetOutputValue( oMesh, id2 );

            if( val1 > val2 )
              {
              OutputPointType temp_pt = q1;
              q1 = q2;
              q2 = temp_pt;
              }

            OutputPixelType temp =
                static_cast< OutputPixelType >(
                  this->Solve( oMesh, iNode, p, id1, q1, id2, q2 ) );

            outputPixel = vnl_math_min( outputPixel, temp );
            }
          }
        else
          {
          // throw one exception here
          }

        qe_it = qe_it2;
        }
      while( qe_it != qe );
      }
    else
      {
      // throw one exception
      }
    }


  OutputVectorRealType
  Solve( OutputMeshType* oDomain,
                const NodeType& iId,
                OutputPointType iCurrentPoint,
                const NodeType& iId1,
                OutputPointType iP1,
                const NodeType& iId2,
                OutputPointType iP2 )
    {
    OutputVectorType Edge1 = iP1 - iCurrentPoint;
    OutputVectorType Edge2 = iP2 - iCurrentPoint;

    OutputVectorRealType sq_norm1 = Edge1.GetSquaredNorm(); //b
    OutputVectorRealType norm1 = vcl_sqrt( sq_norm1 );

    OutputVectorRealType epsilon =
        NumericTraits< OutputVectorRealType >::epsilon();

    if( norm1 > epsilon )
      {
      OutputVectorRealType inv_norm1 = 1. / norm1;
      Edge1 *= inv_norm1;
      }

    OutputVectorRealType sq_norm2 = Edge2.GetSquaredNorm(); //a
    OutputVectorRealType norm2 = vcl_sqrt( sq_norm2 );

    if( norm2 > epsilon )
      {
      OutputVectorRealType inv_norm2 = 1. / norm2;
      Edge2 *= inv_norm2;
      }

    OutputVectorRealType dist1 = iP1.GetVectorFromOrigin().GetNorm();
    OutputVectorRealType dist2 = iP2.GetVectorFromOrigin().GetNorm();

    bool Usable1 = ( this->GetLabelValueForGivenNode( iId1 ) != Superclass::Far );
    bool Usable2 = ( this->GetLabelValueForGivenNode( iId2 ) != Superclass::Far );

    InputPixelType F;
    this->GetInput()->GetPointData( iId, &F );

    if( !Usable1 && Usable2 )
      {
      // only one point is a contributor
      return dist2 + norm2 * F;
      }
    if( Usable1 && !Usable2 )
      {
      // only one point is a contributor
      return dist1 + norm1 * F;
      }

    if( Usable1 && Usable2 )
      {
      OutputVectorRealType dot =
          static_cast< OutputVectorRealType >( Edge1 * Edge2 );

      if( dot > 0. )
        {
        return ComputeUpdate( dist1, dist2,
                             norm2, sq_norm2,
                             norm1, sq_norm1, dot, F );
        }
      else
        {
        // throw an exception here!
        // angle is obtuse, some preprocessing must be done on the input mesh
        }


      }

    return this->m_LargeValue;
    }

  OutputVectorRealType
  ComputeUpdate(
    const OutputVectorRealType& iDist1, const OutputVectorRealType& iDist2,
    const OutputVectorRealType& iNorm1, const OutputVectorRealType& iSqNorm1,
    const OutputVectorRealType& iNorm2, const OutputVectorRealType& iSqNorm2,
    const OutputVectorRealType& iDot, const OutputVectorRealType& iF )
    const
  {
    OutputVectorRealType large_value =
        static_cast< OutputVectorRealType >( this->m_LargeValue );

    OutputVectorRealType t = large_value;

    OutputVectorRealType CosAngle = iDot;
    OutputVectorRealType SinAngle = vcl_sqrt( 1 - iDot * iDot );

    OutputVectorRealType u = iDist2 - iDist1;

    OutputVectorRealType sq_u = u * u;
    OutputVectorRealType f2 = iSqNorm1 + iSqNorm2 - 2. * iNorm1 * iNorm2 * CosAngle;
    OutputVectorRealType f1 = iNorm2 * u * ( iNorm1 * CosAngle - iNorm2 );
    OutputVectorRealType f0 = iSqNorm2 * ( sq_u - iF * iF * iSqNorm1 * SinAngle * SinAngle );

    OutputVectorRealType delta = f1 * f1 - f0 * f2;

    OutputVectorRealType epsilon =
        NumericTraits< OutputVectorRealType >::epsilon();

    if( delta >= 0. )
      {
      if( vnl_math_abs( f2 ) > epsilon )
        {
        t = ( -f1 - vcl_sqrt( delta ) ) / f2;

        // test if we must must choose the other solution
        if( ( t < u ) ||
            ( iNorm2 * ( t - u ) / t < iNorm1 * CosAngle ) ||
            ( iNorm1 / CosAngle < iNorm2 * ( t - u ) / t ) )
          {
          t = ( -f1 + vcl_sqrt( delta ) ) / f2;
          }
        }
      else
        {
        if( vnl_math_abs( f1 ) > epsilon )
          {
          t = -f0 / f1;
          }
        else
          {
          t = - large_value;
          }
        }
      }
    else
      {
      t = - large_value;
      }

    // choose the update from the 2 vertex only if upwind criterion is met
    if( ( u < t ) &&
        ( iNorm1 * CosAngle < iNorm2 * ( t - u ) / t ) &&
        ( iNorm2 * ( t - u ) / t < iNorm1 / CosAngle ) )
      {
      return t + iDist1;
      }
    else
      {
      return vnl_math_min( iNorm2 * iF + iDist1, iNorm1 * iF + iDist2 );
      }
  }



  bool CheckTopology( OutputMeshType* oDomain,
                      const NodeType& iNode )
    {
    (void) oDomain;
    (void) iNode;

    return true;
    }

  void InitializeOutput( OutputMeshType* oDomain )
    {
    this->CopyInputMeshToOutputMesh();

      {
      OutputPointsContainerPointer points = oDomain->GetPoints();
      OutputPointsContainerIterator p_it = points->Begin();
      OutputPointsContainerIterator p_end = points->End();

      while( p_it != p_end )
        {
        oDomain->SetPointData( p_it->Index(), this->m_LargeValue );
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
        NodeType idx = pointsIter->Value().GetNode();
        OutputPixelType outputPixel = pointsIter->Value().GetValue();

        m_Label[idx] = Superclass::Alive;
        oDomain->SetPointData( idx, outputPixel );

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
        NodeType idx = node_it->Value();

        m_Label[idx] = Superclass::Forbidden;
        oDomain->SetPointData( idx, zero );

        ++node_it;
        }
      }
    if ( this->m_TrialPoints )
      {
      NodePairContainerConstIterator pointsIter = this->m_TrialPoints->Begin();
      NodePairContainerConstIterator pointsEnd = this->m_TrialPoints->End();

      while( pointsIter != pointsEnd )
        {
        NodeType idx = pointsIter->Value().GetNode();
        OutputPixelType outputPixel = pointsIter->Value().GetValue();

        m_Label[idx] = Superclass::InitialTrial;

        oDomain->SetPointData( idx, outputPixel );

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
