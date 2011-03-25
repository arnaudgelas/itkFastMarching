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

#ifndef __itkFastMarchingQuadEdgeMeshFilterBase_txx
#define __itkFastMarchingQuadEdgeMeshFilterBase_txx

#include "itkFastMarchingQuadEdgeMeshFilterBase.h"
#include "vnl/vnl_math.h"

namespace itk
{
template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
  TOutputPixel, TOutputMeshTraits >
::FastMarchingQuadEdgeMeshFilterBase() : Superclass()
{}

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
  TOutputPixel, TOutputMeshTraits >
::~FastMarchingQuadEdgeMeshFilterBase()
{}

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
IdentifierType
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
  TOutputPixel, TOutputMeshTraits >
::GetTotalNumberOfNodes() const
{
  return this->GetInput()->GetNumberOfPoints();
}

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
void
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
  TOutputPixel, TOutputMeshTraits >
::SetOutputValue( OutputMeshType* oMesh,
                  const NodeType& iNode,
                  const OutputPixelType& iValue )
{
  oMesh->SetPointData( iNode, iValue );
}

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
const typename
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
TOutputPixel, TOutputMeshTraits >
::OutputPixelType
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
TOutputPixel, TOutputMeshTraits >
::GetOutputValue( OutputMeshType* oMesh, const NodeType& iNode ) const
{
  OutputPixelType outputValue;
  oMesh->GetPointData( iNode, &outputValue );
  return outputValue;
}


template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
const unsigned char
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
TOutputPixel, TOutputMeshTraits >
::GetLabelValueForGivenNode( const NodeType& iNode ) const
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

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
void
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
  TOutputPixel, TOutputMeshTraits >
::SetLabelValueForGivenNode( const NodeType& iNode,
                             const LabelType& iLabel )
{
   m_Label[iNode] = iLabel;
}

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
void
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
  TOutputPixel, TOutputMeshTraits >
::UpdateNeighbors( OutputMeshType* oMesh,
                   const NodeType& iNode )
{
  OutputPointType p;
  oMesh->GetPoint( iNode, &p );

  OutputQEType* qe = p.GetEdge();

  if( qe )
    {
    OutputQEType *qe_it = qe;

    do
      {
      if( qe_it )
        {
        OutputPointIdentifierType neigh_id = qe_it->GetDestination();

        const char label = this->GetLabelValueForGivenNode( neigh_id );

        if ( ( label != Superclass::Alive ) &&
             ( label != Superclass::InitialTrial ) &&
             ( label != Superclass::Forbidden ) )
          {
          this->UpdateValue( oMesh, neigh_id );
          }
        }
      else
        {
        itkGenericExceptionMacro( <<"qe_it is NULL" );
        }
      qe_it = qe_it->GetOnext();
      }
    while( qe_it != qe );
    }
  else
    {
    itkGenericExceptionMacro( <<"qe is NULL" );
    }
}

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
void
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
  TOutputPixel, TOutputMeshTraits >
::UpdateValue( OutputMeshType* oMesh,
               const NodeType& iNode )
  {
  OutputPointType p;
  oMesh->GetPoint( iNode, &p );

  InputPixelType F;
  this->GetInput()->GetPointData( iNode, &F );

  if( F < 0. )
    {
    itkGenericExceptionMacro( << "F < 0." );
    }

  OutputQEType* qe = p.GetEdge();

  OutputPixelType outputPixel = this->m_LargeValue;

  if( qe )
    {
    OutputQEType *qe_it = qe;
    qe_it = qe_it->GetOnext();

    do
      {
      OutputQEType *qe_it2 = qe_it->GetOnext();

      if( qe_it2 )
        {
        if( qe_it->GetLeft() != OutputMeshType::m_NoFace )
          {
          OutputPointIdentifierType id1 = qe_it->GetDestination();
          OutputPointIdentifierType id2 = qe_it2->GetDestination();

          const LabelType label1 =
              static_cast< LabelType >( this->GetLabelValueForGivenNode( id1 ) );
          const LabelType label2 =
              static_cast< LabelType >( this->GetLabelValueForGivenNode( id2 ) );

          bool IsFar1 = ( label1 != Superclass::Far );
          bool IsFar2 = ( label2 != Superclass::Far );

          if( IsFar1 || IsFar2 )
            {
            OutputPointType q1 = oMesh->GetPoint( id1 );
            OutputPointType q2 = oMesh->GetPoint( id2 );

            OutputVectorRealType val1 =
                static_cast< OutputVectorRealType >(
                  this->GetOutputValue( oMesh, id1 ) );

            OutputVectorRealType val2 =
                static_cast< OutputVectorRealType >(
                  this->GetOutputValue( oMesh, id2 ) );

            if( val1 > val2 )
              {
              OutputPointType temp_pt( q1 );
              q1 = q2;
              q2 = temp_pt;

              std::swap( val1, val2 );
              std::swap( IsFar1, IsFar2 );
              std::swap( id1, id2 );
              }

            const OutputVectorRealType temp =
                this->Solve( oMesh, iNode, p, F,
                            id1, q1, IsFar1, val1,
                            id2, q2, IsFar2, val2 );

            outputPixel =
                vnl_math_min( outputPixel,
                              static_cast< OutputPixelType >( temp ) );
            }
          }

        qe_it = qe_it2;
        }
      else
        {
        // throw one exception here
        itkGenericExceptionMacro( << "qe_it2 is NULL" );
        }
      }
    while( qe_it != qe );

    if( outputPixel < this->m_LargeValue )
      {
      this->SetOutputValue( oMesh, iNode, outputPixel );

      this->SetLabelValueForGivenNode( iNode, Superclass::Trial );

      this->m_Heap.push( NodePairType( iNode, outputPixel ) );
      }
    }
  else
    {
    // throw one exception
    itkGenericExceptionMacro( << "qe_it is NULL" );
    }
  }

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
const typename
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
TOutputPixel, TOutputMeshTraits >
::OutputVectorRealType
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
TOutputPixel, TOutputMeshTraits >
::Solve( OutputMeshType* oMesh,
       const NodeType& iId, const OutputPointType& iCurrentPoint,
       const OutputVectorRealType& iF,
       const NodeType& iId1, const OutputPointType& iP1,
       const bool& iIsFar1, const OutputVectorRealType iVal1,
       const NodeType& iId2, const OutputPointType& iP2,
       const bool& iIsFar2, const OutputVectorRealType& iVal2 ) const
  {
  OutputVectorType Edge1 = iP1 - iCurrentPoint;
  OutputVectorType Edge2 = iP2 - iCurrentPoint;

  OutputVectorRealType sq_norm1 = Edge1.GetSquaredNorm();
  OutputVectorRealType norm1 = 0.;

  OutputVectorRealType epsilon =
      NumericTraits< OutputVectorRealType >::epsilon();

  if( sq_norm1 > epsilon )
    {
    norm1 = vcl_sqrt( sq_norm1 );

    OutputVectorRealType inv_norm1 = 1. / norm1;
    Edge1 *= inv_norm1;
    }

  OutputVectorRealType sq_norm2 = Edge2.GetSquaredNorm();
  OutputVectorRealType norm2 = 0.;

  if( sq_norm2 > epsilon )
    {
    norm2 = vcl_sqrt( sq_norm2 );

    OutputVectorRealType inv_norm2 = 1. / norm2;
    Edge2 *= inv_norm2;
    }

  if( !iIsFar1 && iIsFar2 )
    {
    // only one point is a contributor
    return iVal2 + norm2 * iF;
    }
  if( iIsFar1 && !iIsFar2 )
    {
    // only one point is a contributor
    return iVal1 + norm1 * iF;
    }

  if( iIsFar1 && iIsFar2 )
    {
    OutputVectorRealType dot =
        static_cast< OutputVectorRealType >( Edge1 * Edge2 );

    if( dot >= 0. )
      {
      return ComputeUpdate( iVal1, iVal2,
                           norm2, sq_norm2,
                           norm1, sq_norm1, dot, iF );
      }
    else
      {
      OutputVectorRealType sq_norm3, norm3, dot1, dot2;
      OutputPointIdentifierType new_id;

      bool unfolded =
        UnfoldTriangle( oMesh, iId, iCurrentPoint, iId1, iP1,
                        iId2, iP2, norm3, sq_norm3, dot1, dot2 , new_id );

      if( unfolded )
        {
        OutputVectorRealType sq_norm3 = norm3 * norm3;

        OutputVectorRealType val3 =
            static_cast< OutputVectorRealType >(
              this->GetOutputValue( oMesh, new_id ) );
        OutputVectorRealType t1 = ComputeUpdate( iVal1, val3, norm3, sq_norm3,
                                                norm1, sq_norm1, dot1, iF );
        OutputVectorRealType t2 = ComputeUpdate( iVal2, val3, norm3, sq_norm3,
                                                norm2, sq_norm2, dot2, iF );

        return vnl_math_min( t1, t2 );
        }
      else
        {
        return ComputeUpdate( iVal1, iVal2,
                             norm2, sq_norm2,
                             norm1, sq_norm1, dot, iF );
        }
      }


    }

  return this->m_LargeValue;
  }

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
const typename
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
TOutputPixel, TOutputMeshTraits >
::OutputVectorRealType
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
TOutputPixel, TOutputMeshTraits >
::ComputeUpdate(
  const OutputVectorRealType& iVal1, const OutputVectorRealType& iVal2,
  const OutputVectorRealType& iNorm1, const OutputVectorRealType& iSqNorm1,
  const OutputVectorRealType& iNorm2, const OutputVectorRealType& iSqNorm2,
  const OutputVectorRealType& iDot, const OutputVectorRealType& iF )
  const
{
  OutputVectorRealType large_value =
      static_cast< OutputVectorRealType >( this->m_LargeValue );

  OutputVectorRealType t = large_value;

  OutputVectorRealType CosAngle = iDot;
  OutputVectorRealType SinAngle = vcl_sqrt( 1. - iDot * iDot );

  OutputVectorRealType u = iVal2 - iVal1;

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
    return t + iVal1;
    }
  else
    {
    return vnl_math_min( iNorm2 * iF + iVal1, iNorm1 * iF + iVal2 );
    }
}

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
bool
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
TOutputPixel, TOutputMeshTraits >
::UnfoldTriangle(
  OutputMeshType* oMesh,
  const OutputPointIdentifierType& iId, const OutputPointType& iP,
  const OutputPointIdentifierType& iId1, const OutputPointType& iP1,
  const OutputPointIdentifierType& iId2, const OutputPointType &iP2,
  OutputVectorRealType& oNorm, OutputVectorRealType& oSqNorm,
  OutputVectorRealType& oDot1, OutputVectorRealType& oDot2,
  OutputPointIdentifierType& oId ) const
  {
  OutputVectorType Edge1 = iP1 - iP;
  OutputVectorRealType Norm1 = Edge1.GetNorm();
  Edge1 /= Norm1;

  OutputVectorType Edge2 = iP2 - iP;
  OutputVectorRealType Norm2 = Edge2.GetNorm();
  Edge2 /= Norm2;

  OutputVectorRealType dot =
      static_cast< OutputVectorRealType >( Edge1 * Edge2 );

  // the equation of the lines defining the unfolding region [e.g. line 1 : {x ; <x,eq1>=0} ]
  typedef Vector< OutputVectorRealType, 2 > Vector2DType;
  typedef Matrix< OutputVectorRealType, 2, 2 > Matrix2DType;

  Vector2DType v1;
  v1[0] = dot;
  v1[1] = vcl_sqrt( 1. - dot * dot );

  Vector2DType v2;
  v2[0] = 1.;
  v2[1] = 0.;

  Vector2DType x1;
  x1[0] = Norm1;
  x1[1] = 0.;

  Vector2DType x2 = Norm2 * v1;


  // keep track of the starting point
  Vector2DType x_start1( x1 );
  Vector2DType x_start2( x2 );

  OutputPointIdentifierType id1 = iId1;
  OutputPointIdentifierType id2 = iId2;

  OutputQEType *qe = oMesh->FindEdge( id1, id2 );
  qe = qe->GetSym();

  OutputPointType t_pt1 = iP1;
  OutputPointType t_pt2 = iP2;
  OutputPointType t_pt;

  unsigned int nNum = 0;
  while( nNum<50 && qe->GetLeft() != OutputMeshType::m_NoFace )
    {
    OutputQEType* qe_Lnext = qe->GetLnext();
    OutputPointIdentifierType t_id = qe_Lnext->GetDestination();

    oMesh->GetPoint( t_id, &t_pt );

    Edge1 = t_pt2 - t_pt1;
    Norm1 = Edge1.GetNorm();
    Edge1 /= Norm1;

    Edge2 = t_pt - t_pt1;
    Norm2 = Edge2.GetNorm();
    Edge2 /= Norm2;

    /* compute the position of the new point x on the unfolding plane (via a rotation of -alpha on (x2-x1)/rNorm1 )
            | cos(alpha) sin(alpha)|
        x = |-sin(alpha) cos(alpha)| * [x2-x1]*rNorm2/rNorm1 + x1   where cos(alpha)=dot
    */
    Vector2DType vv = (x2 - x1) * Norm2 / Norm1;

    dot = Edge1 * Edge2;

    Matrix2DType rotation;
    rotation[0][0] = dot;
    rotation[0][1] = vcl_sqrt( 1. - dot * dot );
    rotation[1][0] = - rotation[0][1];
    rotation[1][1] = dot;

    Vector2DType x = rotation * vv + x1;


    /* compute the intersection points.
       We look for x=x1+lambda*(x-x1) or x=x2+lambda*(x-x2) with <x,eqi>=0, so */
    OutputVectorRealType lambda11 = - ( x1 * v1 ) / ( ( x - x1 ) * v1 );	// left most
    OutputVectorRealType lambda12 = - ( x1 * v2 ) / ( ( x - x1 ) * v2 );	// right most
    OutputVectorRealType lambda21 = - ( x2 * v1 ) / ( ( x - x2 ) * v1 );	// left most
    OutputVectorRealType lambda22 = - ( x2 * v2 ) / ( ( x - x2 ) * v2 );	// right most
    bool bIntersect11 = (lambda11>=0.) && (lambda11<=1.);
    bool bIntersect12 = (lambda12>=0.) && (lambda12<=1.);
    bool bIntersect21 = (lambda21>=0.) && (lambda21<=1.);
    bool bIntersect22 = (lambda22>=0.) && (lambda22<=1.);

    if( bIntersect11 && bIntersect12 )
      {
      qe = oMesh->FindEdge( id1, t_id );
      qe->GetSym();

      t_pt2 = t_pt;
      x2 = x;
      }
    else
      {
      if( bIntersect21 && bIntersect22 )
        {
        qe = oMesh->FindEdge( id2, t_id );
        qe->GetSym();

        t_pt1 = t_pt;
        x1 = x;
        }
      else
        { // bIntersect11 && !bIntersect12 && !bIntersect21 && bIntersect22

        // that's it, we have found the point
        oSqNorm = x.GetSquaredNorm();

        if( oSqNorm > NumericTraits< OutputVectorRealType >::epsilon() )
          {
          oNorm = vcl_sqrt( oSqNorm );
          oDot1 = x * x_start1 / ( oNorm * x_start1.GetNorm() );
          oDot2 = x * x_start2 / ( oNorm * x_start2.GetNorm() );
          }
        else
          {
          oNorm = 0.;
          oDot1 = 0.;
          oDot2 = 0.;
          }

        oId = t_id;
        return true;
        }
      }
    ++nNum;
    }

    return false;
  }

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
bool
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
TOutputPixel, TOutputMeshTraits >
::CheckTopology( OutputMeshType* oMesh,
                    const NodeType& iNode )
  {
  (void) oMesh;
  (void) iNode;

  return true;
  }

template< unsigned int VDimension,
          typename TInputPixel,
          class TInputMeshTraits,
          typename TOutputPixel,
          class TOutputMeshTraits >
void
FastMarchingQuadEdgeMeshFilterBase< VDimension, TInputPixel, TInputMeshTraits,
TOutputPixel, TOutputMeshTraits >
::InitializeOutput( OutputMeshType* oMesh )
  {
  this->CopyInputMeshToOutputMeshGeometry();

  // Check that the input mesh is made of triangles
    {
    OutputCellsContainerPointer cells = oMesh->GetCells();
    OutputCellsContainerConstIterator c_it = cells->Begin();
    OutputCellsContainerConstIterator c_end = cells->End();

    while( c_it != c_end )
      {
      OutputCellType* cell = c_it.Value();

      if( cell->GetNumberOfPoints() != 3 )
        {
        itkGenericExceptionMacro( << "Input mesh has non triangular faces" );
        }
      ++c_it;
      }
    }

  OutputPointsContainerPointer points = oMesh->GetPoints();

  OutputPointDataContainerPointer pointdata =
      OutputPointDataContainer::New();
  pointdata->Reserve( points->Size() );

  OutputPointsContainerIterator p_it = points->Begin();
  OutputPointsContainerIterator p_end = points->End();

  while( p_it != p_end )
    {
    pointdata->SetElement( p_it->Index(), this->m_LargeValue );
    ++p_it;
    }

  oMesh->SetPointData( pointdata );

  m_Label.clear();

  if ( this->m_AlivePoints )
    {
    NodePairContainerConstIterator pointsIter = this->m_AlivePoints->Begin();
    NodePairContainerConstIterator pointsEnd = this->m_AlivePoints->End();

    while( pointsIter != pointsEnd )
      {
      // get node from alive points container
      NodeType idx = pointsIter->Value().GetNode();

      if( points->IndexExists( idx ) )
        {
        OutputPixelType outputPixel = pointsIter->Value().GetValue();

        this->SetLabelValueForGivenNode( idx, Superclass::Alive );
        this->SetOutputValue( oMesh, idx, outputPixel );
        }

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

      if( points->IndexExists( idx ) )
        {
        this->SetLabelValueForGivenNode( idx, Superclass::Forbidden );
        this->SetOutputValue( oMesh, idx, zero );
        }

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

      if( points->IndexExists( idx ) )
        {
        OutputPixelType outputPixel = pointsIter->Value().GetValue();

        this->SetLabelValueForGivenNode( idx, Superclass::InitialTrial );
        this->SetOutputValue( oMesh, idx, outputPixel );

        this->m_Heap.push( pointsIter->Value() );
        }

      ++pointsIter;
      }
    }
  }
}

#endif // __itkFastMarchingQuadEdgeMeshFilterBase_txx
