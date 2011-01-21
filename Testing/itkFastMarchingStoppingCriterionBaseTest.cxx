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

#include "itkFastMarchingStoppingCriterionBase.h"
#include "itkFastMarchingBase.h"

namespace itk
{
  template< class TTraits >
  class FastMarchingStoppingCriterionBaseHelperTest :
      public FastMarchingStoppingCriterionBase< TTraits >
    {
  public:
    typedef FastMarchingStoppingCriterionBaseHelperTest Self;
    typedef FastMarchingStoppingCriterionBase< TTraits > Superclass;
    typedef SmartPointer< Self >              Pointer;
    typedef SmartPointer< const Self >        ConstPointer;

    typedef typename Superclass::NodeType NodeType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(FastMarchingStoppingCriterionBaseHelperTest,
                 FastMarchingStoppingCriterionBase );

    void SetCurrentNode( const NodeType& iNode ) { (void) iNode; }
    bool IsSatisfied() const { return true; }

  protected:
    FastMarchingStoppingCriterionBaseHelperTest() : Superclass() {}
    ~FastMarchingStoppingCriterionBaseHelperTest() {}
  private:
    FastMarchingStoppingCriterionBaseHelperTest( const Self& );
    void operator = ( const Self& );
    };
}

int main( int argc, char* argv[] )
  {
  (void) argc;
  (void) argv;

   typedef itk::ImageFastMarchingTraits< 2, float, float >
    ImageTraits;
  typedef itk::FastMarchingStoppingCriterionBaseHelperTest< ImageTraits >
    ImageStoppingCriterionType;

  ImageStoppingCriterionType::Pointer image_criterion =
      ImageStoppingCriterionType::New();

  typedef itk::MeshFastMarchingTraits< 3,
      float,
      itk::DefaultStaticMeshTraits< float, 3, 3 >,
      float,
      itk::DefaultStaticMeshTraits< float, 3, 3 > >
    MeshTraits;

  typedef itk::FastMarchingStoppingCriterionBaseHelperTest< MeshTraits >
    MeshStoppingCriterionType;

  MeshStoppingCriterionType::Pointer mesh_criterion =
      MeshStoppingCriterionType::New();

  return EXIT_SUCCESS;
  }
