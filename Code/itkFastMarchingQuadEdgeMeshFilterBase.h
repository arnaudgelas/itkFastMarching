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
  typedef typename OutputMeshType::PointDataContainer
                                                    OutputPointDataContainer;
  typedef typename OutputPointDataContainer::Pointer
                                                    OutputPointDataContainerPointer;

  typedef typename OutputMeshType::CellsContainer   OutputCellsContainer;
  typedef typename OutputCellsContainer::Pointer    OutputCellsContainerPointer;
  typedef typename OutputCellsContainer::ConstIterator
                                                    OutputCellsContainerConstIterator;
  typedef typename OutputMeshType::CellType         OutputCellType;


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

  void SetOutputValue( OutputMeshType* oDomain,
                      const NodeType& iNode,
                      const OutputPixelType& iValue )
    {
    oDomain->SetPointData( iNode, iValue );
    }

  const OutputPixelType GetOutputValue( OutputMeshType* oDomain,
                                  const NodeType& iNode ) const
    {
    OutputPixelType outputValue;
    oDomain->GetPointData( iNode, &outputValue );
    return outputValue;
    }

  const unsigned char GetLabelValueForGivenNode( const NodeType& iNode ) const
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



  void UpdateValue( OutputMeshType* oMesh,
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


  const OutputVectorRealType
  Solve( OutputMeshType* oDomain,
         const NodeType& iId, const OutputPointType& iCurrentPoint,
         const OutputVectorRealType& iF,
         const NodeType& iId1, const OutputPointType& iP1,
         const bool& iIsFar1, const OutputVectorRealType iVal1,
         const NodeType& iId2, const OutputPointType& iP2,
         const bool& iIsFar2, const OutputVectorRealType& iVal2 )
  const
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
        // throw an exception here!
        // angle is obtuse, some preprocessing must be done on the input mesh
        itkWarningMacro(
              <<"Not yet implemented for meshes with obtuse angle." << std::endl
              << iId <<" iCurrentPoint: " << iCurrentPoint <<std::endl
              << iId1 <<" iP1: " << iP1 << std::endl
              << iId2 <<" iP2: " << iP2 << std::endl
              <<"dot = " << dot << std::endl
              );
        }


      }

    return this->m_LargeValue;
    }

  const OutputVectorRealType
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
    OutputVectorRealType SinAngle = vcl_sqrt( 1. - iDot * iDot );

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
    this->CopyInputMeshToOutputMeshGeometry();

    // Check that the input mesh is made of triangles
      {
      OutputCellsContainerPointer cells = oDomain->GetCells();
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

    OutputPointsContainerPointer points = oDomain->GetPoints();

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

    oDomain->SetPointData( pointdata );

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
          this->SetOutputValue( oDomain, idx, outputPixel );
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
          this->SetOutputValue( oDomain, idx, zero );
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
          this->SetOutputValue( oDomain, idx, outputPixel );

          this->m_Heap.push( pointsIter->Value() );
          }

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
