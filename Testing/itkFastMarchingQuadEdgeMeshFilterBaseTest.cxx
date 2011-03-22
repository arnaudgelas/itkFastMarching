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

  typedef itk::QuadEdgeMesh< CoordType, Dimension, Traits > MeshType;

  typedef itk::FastMarchingQuadEdgeMeshFilterBase< Dimension, PixelType, Traits,
      PixelType, Traits > FastMarchingType;
  FastMarchingType::Pointer filter = FastMarchingType::New();

  return EXIT_SUCCESS;
}
