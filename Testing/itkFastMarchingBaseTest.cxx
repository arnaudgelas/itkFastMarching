#include "itkFastMarchingBase.h"
#include "itkDefaultStaticMeshTraits.h"

namespace itk
{
template< class TTraits >
class FastMarchingBaseTestHelper : public FastMarchingBase< TTraits >
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

  typedef typename Superclass::NodeContainerType  NodeContainerType;
  typedef typename Superclass::NodeType           NodeType;

  typedef typename Superclass::OutputPixelType  OutputPixelType;
  typedef typename Superclass::LabelType        LabelType;

  void SetAliveNodes( NodeContainerType iNodes )
    { (void) iNodes; }
  virtual void AddAliveNode( NodeType iNode, OutputPixelType iValue )
    { (void) iNode; (void) iValue; }

  virtual void SetForbiddenNodes( NodeType iNodes )
    { (void) iNodes; }
  virtual void AddForbiddenNode( NodeType iNode )
    { (void) iNode; }

protected:
  FastMarchingBaseTestHelper() {}
  ~FastMarchingBaseTestHelper() {}

  LabelType GetLabelValueForGivenNode( NodeType iNode )
    {
    return Superclass::Far;
    }

  void SetLabelValueForGivenNode( NodeType iNode, LabelType iLabel )
    { (void) iNode; (void) iLabel; }

  void UpdateNeighbors( NodeType iNode )
    { (void) iNode; }

  void CheckTopology( NodeType iNode )
    { (void) iNode; }

private:
  FastMarchingBaseTestHelper( const Self& );
  void operator = ( const Self& );
};
}

int main( int argc, char* argv[] )
{
  if( argc != 2 )
    {
    return EXIT_FAILURE;
    }

  if( atoi( argv[1] ) == 0 )
    {
    typedef itk::ImageFastMarchingTraits< 2, float, float >
      ImageTraits;
    typedef ImageTraits::InputDomainType InputImageType;
    InputImageType::Pointer input = InputImageType::New();

    typedef itk::FastMarchingBaseTestHelper< ImageTraits > ImageFastMarching;
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
      typedef itk::MeshFastMarchingTraits< 3,
          float,
          itk::DefaultStaticMeshTraits< float, 3, 3 >,
          float,
          itk::DefaultStaticMeshTraits< float, 3, 3 > >
        MeshTraits;
      typedef MeshTraits::InputDomainType InputMeshType;
      InputMeshType::Pointer input = InputMeshType::New();

      typedef itk::FastMarchingBaseTestHelper< MeshTraits > MeshFastMarching;
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
