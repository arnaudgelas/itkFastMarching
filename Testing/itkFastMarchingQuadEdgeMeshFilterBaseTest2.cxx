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

#include "itkFastMarchingQuadEdgeMeshFilterBase.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshExtendedTraits.h"
#include "itkRegularSphereMeshSource.h"
#include "itkFastMarchingThresholdStoppingCriterion.h"
#include "itkQuadEdgeMeshScalarDataVTKPolyDataWriter.h"

int main( int argc, char* argv[] )
{
  typedef float PixelType;
  typedef double CoordType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMeshExtendedTraits <
    PixelType,  // type of data for vertices
    Dimension,  // geometrical dimension of space
    2,          // Mac topological dimension of a cell
    CoordType,  // type for point coordinate
    CoordType,  // type for interpolation weight
    PixelType,  // type of data for cell
    bool,       // type of data for primal edges
    bool        // type of data for dual edges
  > Traits;

  typedef itk::FastMarchingQuadEdgeMeshFilterBase< Dimension, PixelType, Traits,
      PixelType, Traits > FastMarchingType;
  typedef FastMarchingType::InputMeshType MeshType;

  typedef MeshType::PointsContainer PointsContainer;
  typedef PointsContainer::Pointer  PointsContainerPointer;

  typedef MeshType::PointDataContainer PointDataContainer;
  typedef PointDataContainer::Pointer  PointDataContainerPointer;

  typedef MeshType::CellType CellType;
  typedef itk::QuadEdgeMeshPolygonCell< CellType > QEPolygonCellType;


  // Let's create here a plane!
  MeshType::Pointer plane = MeshType::New();

  PointsContainerPointer points = PointsContainer::New();
  PointDataContainerPointer pointdata = PointDataContainer::New();

  points->Reserve( 100 );
  pointdata->Reserve( 100 );

  MeshType::PointType p;
  p[2] = 0.;

  int k = 0;

  for( int i = 0; i < 10; i++ )
    {
    p[0] = static_cast< CoordType >( i );

    for( int j = 0; j < 10; j++ )
      {
      p[1] = static_cast< CoordType >( j );
      points->SetElement( k, p );
      pointdata->SetElement( k, 1. );
      k++;
      }
    }

  plane->SetPoints( points );
  plane->SetPointData( pointdata );

  k = 0;

  for( int i = 0; i < 9; i++ )
    {
    for( int j = 0; j < 9; j++ )
      {
      plane->AddFaceTriangle( k, k+1, k+11 );
      plane->AddFaceTriangle( k, k+11, k+10 );
      k++;
      }
    k++;
    }

  typedef FastMarchingType::NodeType NodeType;
  typedef FastMarchingType::NodePairType NodePairType;
  typedef FastMarchingType::NodeContainerType NodeContainerType;
  typedef FastMarchingType::NodePairContainerType NodePairContainerType;

  NodePairContainerType::Pointer trial = NodePairContainerType::New();

  NodePairType node_pair( 0, 0. );
  trial->push_back( node_pair );

  typedef itk::FastMarchingThresholdStoppingCriterion< NodeType, PixelType >
      CriterionType;
  CriterionType::Pointer criterion = CriterionType::New();
  criterion->SetThreshold( 100. );

  FastMarchingType::Pointer fmm_filter = FastMarchingType::New();
  fmm_filter->SetInput( plane );
  fmm_filter->SetTrialPoints( trial );
  fmm_filter->SetStoppingCriterion( criterion );

  try
    {
    fmm_filter->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::QuadEdgeMeshScalarDataVTKPolyDataWriter< MeshType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( fmm_filter->GetOutput() );
  writer->SetFileName( "output.vtk" );
  writer->Update();


  return EXIT_SUCCESS;
}
