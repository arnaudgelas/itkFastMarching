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

#include "itkFastMarchingTraits.h"
#include "itkFastMarchingBase.h"
#include "itkDefaultStaticMeshTraits.h"

namespace itk
{
template< class TTraits >
class FastMarchingBaseTestHelper :
    public FastMarchingBase< TTraits >
{
public:
  typedef FastMarchingBaseTestHelper Self;
  typedef FastMarchingBase< TTraits > Superclass;
  typedef SmartPointer< Self > Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FastMarchingBaseTestHelper, FastMarchingBase);

  typedef typename Superclass::OutputDomainType OutputDomainType;

  typedef typename Superclass::NodeContainerType  NodeContainerType;
  typedef typename Superclass::NodeType           NodeType;

  typedef typename Superclass::OutputPixelType  OutputPixelType;
  typedef typename Superclass::LabelType        LabelType;

protected:
  FastMarchingBaseTestHelper() {}
  ~FastMarchingBaseTestHelper() {}

  IdentifierType GetTotalNumberOfNodes() const
    { return 1; }

  OutputPixelType GetOutputValue( OutputDomainType* ,
                                  const NodeType& ) const
    {
    return NumericTraits< OutputPixelType >::Zero;
    }

  char GetLabelValueForGivenNode( const NodeType& ) const
    {
    return Superclass::Far;
    }

  void SetLabelValueForGivenNode( const NodeType& ,
                                 const LabelType& )
    {}

  void UpdateNeighbors( OutputDomainType* , const NodeType& )
    {}

  void UpdateValue( OutputDomainType* , const NodeType& )
    {}

  bool CheckTopology( OutputDomainType* , const NodeType&  )
    { return true; }

  void InitializeOutput( OutputDomainType* ) {}

private:
  FastMarchingBaseTestHelper( const Self& );
  void operator = ( const Self& );
};
}






// -----------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
  if( argc != 2 )
    {
    return EXIT_FAILURE;
    }

  const unsigned Dimension = 3;
  typedef float PixelType;

  bool exception_caught = false;

  if( atoi( argv[1] ) == 0 )
    {
    typedef itk::ImageFastMarchingTraits< Dimension, PixelType, PixelType >
      ImageTraits;
    typedef ImageTraits::InputDomainType InputImageType;
    typedef ImageTraits::NodeType ImageNodeType;

    typedef itk::FastMarchingStoppingCriterionBase< ImageNodeType, PixelType >
        ImageCriterionType;

    InputImageType::Pointer input = InputImageType::New();

    typedef itk::FastMarchingBaseTestHelper< ImageTraits >
        ImageFastMarching;
    ImageFastMarching::Pointer fmm = ImageFastMarching::New();
    fmm->SetInput( input );

    try
      {
      fmm->Update();
      }
    catch( itk::ExceptionObject & excep )
      {
      std::cerr << "Exception caught !" << std::endl;
      std::cerr << excep << std::endl;
      exception_caught = true;
      }

    typedef ImageFastMarching::OutputDomainType OutputImageType;
    OutputImageType::Pointer output = fmm->GetOutput();

    (void) output;
    }
  else
    {
    if( atoi( argv[1] ) == 1 )
      {
      typedef itk::MeshFastMarchingTraits< Dimension,
          PixelType,
          itk::QuadEdgeMeshTraits< PixelType, 3, bool, bool >,
          PixelType,
          itk::QuadEdgeMeshTraits< PixelType, 3, bool, bool > >
        MeshTraits;
      typedef MeshTraits::NodeType MeshNodetype;

      typedef itk::FastMarchingStoppingCriterionBase< MeshNodetype, PixelType >
          MeshCriterionType;

      typedef MeshTraits::InputDomainType InputMeshType;
      InputMeshType::Pointer input = InputMeshType::New();

      typedef itk::FastMarchingBaseTestHelper< MeshTraits >
          MeshFastMarching;
      MeshFastMarching::Pointer fmm = MeshFastMarching::New();
      fmm->SetInput( input );

      try
        {
        fmm->Update();
        }
      catch( itk::ExceptionObject & excep )
        {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        exception_caught = true;
        }

      typedef MeshFastMarching::OutputDomainType OutputMeshType;
      OutputMeshType::Pointer output = fmm->GetOutput();

      (void) output;
      }
    }

  if( exception_caught )
    {
    return EXIT_SUCCESS;
    }
  else
    {
    return EXIT_FAILURE;
    }
}
