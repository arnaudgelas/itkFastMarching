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

#ifndef __Test_SpeedToPath_GradientDescent_ND
#define __Test_SpeedToPath_GradientDescent_ND

#include "ReadPathFile.h"

/////////////////////////////////////////////////////////////
// Template for SpeedToPath with GradientDescentOptimizer
template <int VDimension>
int Test_SpeedToPath_GradientDescent_ND(int argc, char* argv[])
{
  const unsigned int Dimension = VDimension;
  typedef float PixelType;
  typedef unsigned char OutputPixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  typedef itk::PolyLineParametricPath< Dimension > PathType;
  typedef itk::SpeedFunctionToPathFilter< ImageType, PathType > PathFilterType;
  typedef typename PathFilterType::CostFunctionType::CoordRepType CoordRepType;
  typedef itk::PathIterator< OutputImageType, PathType > PathIteratorType;

  try
    {
    // Print header info
    if ( argc != 6 )
      {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0];
      std::cerr << " OutputFilename";
      std::cerr << " SpeedFilename";
      std::cerr << " PathFilename";
      std::cerr << " TerminationValue";    // Good default = 2.0
      std::cerr << " NumberOfIterations";  // Good default = 1000
      std::cerr << std::endl;
      return EXIT_FAILURE;
      }

    // Get arguments
    unsigned int argi = 1;
    char* OutputFilename = argv[argi++];
    char* SpeedFilename = argv[argi++];
    char* PathFilename = argv[argi++];
    float TerminationValue = atof( argv[argi++] );
    unsigned int NumberOfIterations = atoi( argv[argi++] );
    // NOTE: Points will be read from the command line later

    // Read speed function
    std::cout << "Speed: " << SpeedFilename << std::endl;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( SpeedFilename );
    reader->Update();
    typename ImageType::Pointer speed = reader->GetOutput();
    speed->DisconnectPipeline();

    // Create Interpolator
    typedef itk::LinearInterpolateImageFunction<ImageType, CoordRepType>
        InterpolatorType;
    typename InterpolatorType::Pointer interp = InterpolatorType::New();

    // Create Cost Function
    typename PathFilterType::CostFunctionType::Pointer
        cost = PathFilterType::CostFunctionType::New();
    cost->SetInterpolator( interp );

    // Create GradientDescentOptimizer
    typedef itk::GradientDescentOptimizer OptimizerType;
    typename OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetNumberOfIterations( NumberOfIterations );

    // Create path filter
    typename PathFilterType::Pointer pathFilter = PathFilterType::New();
    pathFilter->SetInput( speed );
    pathFilter->SetCostFunction( cost );
    pathFilter->SetOptimizer( optimizer );
    pathFilter->SetTerminationValue( TerminationValue );

    // Read path file
    if ( ReadPathFile<PathFilterType>(PathFilename, pathFilter) == EXIT_FAILURE )
      {
      std::cerr << "Failed to read path file: " << PathFilename << std::endl;
      return EXIT_FAILURE;
      }

    // Compute the path
    std::cout << "Computing path..." << std::endl;
    itk::TimeProbe time;
    try
      {
      time.Start( );
      pathFilter->Update( );
      time.Stop( );
      }
    catch (itk::ExceptionObject & err)
      {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }

    std::cout << std::setprecision(3) << "Path computed in: ";
    std::cout << time.GetMeanTime() << " seconds" << std::endl;

    // Allocate output image
    typename OutputImageType::Pointer output = OutputImageType::New();
    output->SetRegions( speed->GetLargestPossibleRegion() );
    output->SetSpacing( speed->GetSpacing() );
    output->SetOrigin( speed->GetOrigin() );
    output->Allocate( );
    output->FillBuffer( itk::NumericTraits<OutputPixelType>::Zero );

    // Rasterize path
    std::cout << pathFilter->GetNumberOfOutputs() <<std::endl;

    for (unsigned int i=0; i<pathFilter->GetNumberOfOutputs(); i++)
      {
      std::cout <<"** " << i <<std::endl;
      // Get the path
      typename PathType::Pointer path = pathFilter->GetOutput( i );

      // Check path is valid
      if ( path->GetVertexList()->Size() == 0 )
        {
        std::cout << "WARNING: Path " << (i+1) << " contains no points!" << std::endl;
        continue;
        }
      // Iterate path and convert to image
      std::cout << "Rasterizing path..." << std::endl;
      PathIteratorType it( output, path );
      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
        it.Set( itk::NumericTraits<OutputPixelType>::max() );
        }
      }

    // Write output
    std::cout << "Output: " << OutputFilename << std::endl;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( OutputFilename );
    writer->SetInput( output );
    try
      {
      writer->Update();
      }
    catch (itk::ExceptionObject & err)
      {
      std::cerr << "writer" <<std::endl;
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }
    }
  catch (itk::ExceptionObject & err)
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  //Return
  return EXIT_SUCCESS;
}

#endif // TEST_SPEEDTOPATH_GRADIENTDESCENT_ND_H
