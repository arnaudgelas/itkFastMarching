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

#include "itkFastMarchingBase.h"
#include "itkDefaultStaticMeshTraits.h"

namespace itk
{
template< class TTraits, class TCriterion >
class FastMarchingBaseTestHelper :
    public FastMarchingBase< TTraits, TCriterion >
{
public:
  typedef FastMarchingBaseTestHelper Self;
  typedef FastMarchingBase< TTraits, TCriterion > Superclass;
  typedef SmartPointer< Self > Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FastMarchingBaseTestHelper, FastMarchingBase);

  typedef typename Superclass::NodeContainerType  NodeContainerType;
  typedef typename Superclass::NodeType           NodeType;

  typedef typename Superclass::OutputPixelType  OutputPixelType;
  typedef typename Superclass::LabelType        LabelType;

protected:
  FastMarchingBaseTestHelper() {}
  ~FastMarchingBaseTestHelper() {}

  char GetLabelValueForGivenNode( NodeType iNode )
    {
    return Superclass::Far;
    }

  void SetLabelValueForGivenNode( NodeType iNode, LabelType iLabel )
    { (void) iNode; (void) iLabel; }

  void UpdateNeighbors( NodeType iNode )
    { (void) iNode; }

  void UpdateValue( NodeType iNode )
    { (void) iNode; }

  bool CheckTopology( NodeType iNode )
    { return true; }

  virtual void InitializeOutput() {}

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

  if( atoi( argv[1] ) == 0 )
    {
    typedef itk::ImageFastMarchingTraits< Dimension, PixelType, PixelType >
      ImageTraits;
    typedef ImageTraits::InputDomainType InputImageType;
    typedef ImageTraits::NodeType ImageNodeType;

    typedef itk::FastMarchingStoppingCriterionBase< ImageNodeType, PixelType >
        ImageCriterionType;

    InputImageType::Pointer input = InputImageType::New();

    typedef itk::FastMarchingBaseTestHelper< ImageTraits, ImageCriterionType >
        ImageFastMarching;
    ImageFastMarching::Pointer fmm = ImageFastMarching::New();
    fmm->SetInput( input );
    fmm->Update();

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
          itk::DefaultStaticMeshTraits< PixelType, 3, 3 >,
          PixelType,
          itk::DefaultStaticMeshTraits< PixelType, 3, 3 > >
        MeshTraits;
      typedef MeshTraits::NodeType MeshNodetype;

      typedef itk::FastMarchingStoppingCriterionBase< MeshNodetype, PixelType >
          MeshCriterionType;

      typedef MeshTraits::InputDomainType InputMeshType;
      InputMeshType::Pointer input = InputMeshType::New();

      typedef itk::FastMarchingBaseTestHelper< MeshTraits, MeshCriterionType >
          MeshFastMarching;
      MeshFastMarching::Pointer fmm = MeshFastMarching::New();
      fmm->SetInput( input );
      fmm->Update();

      typedef MeshFastMarching::OutputDomainType OutputMeshType;
      OutputMeshType::Pointer output = fmm->GetOutput();

      (void) output;
      }
    }

  return EXIT_SUCCESS;
}
