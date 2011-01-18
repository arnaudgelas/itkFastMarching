#include "itkFastMarchingBase.h"

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

  void SetTrialNodes( NodeContainerType iNodes )
    { (void) iNodes; }
  void AddTrialNode( NodeType iNode, OutputPixelType iValue )
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
  typedef itk::ImageFastMarchingTraits< 2, float, float >
    ImageTraits;
  typedef itk::FastMarchingBaseTestHelper< ImageTraits > ImageFastMarching;
  ImageFastMarching::Pointer fmm = ImageFastMarching::New();

  return EXIT_SUCCESS;
}
