#include "itkFastMarchingBase.h"

int main( int argc, char* argv[] )
{
  typedef ImageFastMarchingTraits< 2, float, float >
    ImageTraits;
  typedef FastMarchingBase< ImageTraits > ImageFastMarching;
  ImageFastMarching::Pointer fmm = ImageFastMarching::New();

  return EXIT_SUCCESS;
}
