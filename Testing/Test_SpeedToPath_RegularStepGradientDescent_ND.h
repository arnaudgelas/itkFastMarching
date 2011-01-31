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

#ifndef TEST_SPEEDTOPATH_REGULARSTEPGRADIENTDESCENT_ND_H
#define TEST_SPEEDTOPATH_REGULARSTEPGRADIENTDESCENT_ND_H


#include "ReadPathFile.h"

/////////////////////////////////////////////////////////////
// Template for SpeedToPath with RegularStepGradientDescentOptimizer
template <int VDimension>
int Test_SpeedToPath_RegularStepGradientDescent_ND(int argc, char* argv[])
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
    if ( argc != 8 )
      {
      std::cerr << "Usage: " << std::endl;
      std::cerr << argv[0];
      std::cerr << " OutputFilename";
      std::cerr << " SpeedFilename";
      std::cerr << " PathFilename";
      std::cerr << " TerminationValue";    // Good default = 2.0
      std::cerr << " NumberOfIterations";  // Good default = 1000
      std::cerr << " StepLengthFactor";    // Good default = 1.0
      std::cerr << " StepLengthRelax";     // Good default = 0.999
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
    float StepLengthFactor = atof( argv[argi++] );
    float StepLengthRelax = atof( argv[argi++] );

    // NOTE: Points will be read from the command line later
    // Read speed function
    std::cout << "Speed: " << SpeedFilename << std::endl;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( SpeedFilename );
    reader->Update();
    typename ImageType::Pointer speed = reader->GetOutput();
    speed->DisconnectPipeline();

    // Compute the minimum spacing
    typename ImageType::SpacingType spacing = speed->GetSpacing();
    double minspacing = spacing[0];

    for (unsigned int dim=0; dim<Dimension; dim++)
      {
      if (spacing[dim] < minspacing)
        {
        minspacing = spacing[dim];
        }
      }
    // Create Interpolator
    typedef itk::LinearInterpolateImageFunction<ImageType, CoordRepType>
        InterpolatorType;
    typename InterpolatorType::Pointer interp = InterpolatorType::New();

    // Create Cost Function
    typename PathFilterType::CostFunctionType::Pointer cost =
        PathFilterType::CostFunctionType::New();
    cost->SetInterpolator( interp );

    // Create RegularStepGradientDescentOptimizer
    typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
    typename OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetNumberOfIterations( NumberOfIterations );
    optimizer->SetMaximumStepLength( 1.0*StepLengthFactor*minspacing );
    optimizer->SetMinimumStepLength( 0.5*StepLengthFactor*minspacing );
    optimizer->SetRelaxationFactor( StepLengthRelax );

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
    time.Start( );
    pathFilter->Update( );
    time.Stop( );
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
    for (unsigned int i=0; i<pathFilter->GetNumberOfOutputs(); i++)
      {
      // Get the path
      typename PathType::Pointer path = pathFilter->GetOutput( i );

      // Check path is valid
      if ( path->GetVertexList()->Size() == 0 )
        {
        std::cout << "WARNING: Path " << (i+1) << " contains no points!" << std::endl;
        }
      else
        {
        // Iterate path and convert to image
        std::cout << "Rasterizing path..." << std::endl;
        PathIteratorType it( output, path );
        it.GoToBegin();

        while ( !it.IsAtEnd() )
          {
          it.Set( itk::NumericTraits<OutputPixelType>::max() );
          ++it;
          }
        }
      }

    // Write output
    std::cout << "Output: " << OutputFilename << std::endl;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( OutputFilename );
    writer->SetInput( output );
    writer->Update();
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

#endif // TEST_SPEEDTOPATH_REGULARSTEPGRADIENTDESCENT_ND_H
