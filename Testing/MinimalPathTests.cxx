/*=========================================================================

  Filter:    MinimalPath
  Program:   Insight Segmentation & Registration Toolkit
  Module:    MinimalPathTests.cxx
  Language:  C++
  Date:      2008/03/01
  Version:   2.0

  Portions of this code are covered under the ITK and VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
=========================================================================*/

#if defined(_MSC_VER)
//Warning about: identifier was truncated to '255' characters in the debug information (MVC6.0 Debug)
#pragma warning( disable : 4786 )
#endif

// General includes
#include <string>
#include <iostream>
#include <iomanip>
#include "itksys/SystemTools.hxx"

// ITK includes
#include "itkNumericTraits.h"
#include "itkTimeProbe.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPolyLineParametricPath.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkArrivalFunctionToPathFilter.h"
#include "itkSpeedFunctionToPathFilter.h"
#include "itkPathIterator.h"
#include "itkGradientDescentOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkIterateNeighborhoodOptimizer.h"
#include "itkTubeSpatialObject.h"
#include "itkTubeSpatialObjectPoint.h"
#include "itkSpatialObjectPoint.h"
#include "itkSpatialObjectWriter.h"

#include "Test_SpeedToPath_GradientDescent_ND.h"
#include "Test_SpeedToPath_RegularStepGradientDescent_ND.h"
#include "Test_SpeedToPath_IterateNeighborhood_ND.h"



/////////////////////////////////////////////////////////////
// Test for SpeedToPath with GradientDescentOptimizer for 2D images
int Test_SpeedToPath_GradientDescent_2D(int argc, char* argv[])
{
  return Test_SpeedToPath_GradientDescent_ND<2>(argc, argv);
}

/////////////////////////////////////////////////////////////
// Test for SpeedToPath with GradientDescentOptimizer for 3D images
int Test_SpeedToPath_GradientDescent_3D(int argc, char* argv[])
{
  return Test_SpeedToPath_GradientDescent_ND<3>(argc, argv);
}

/////////////////////////////////////////////////////////////
// Test for SpeedToPath with RegularStepGradientDescentOptimizer for 2D images
int Test_SpeedToPath_RegularStepGradientDescent_2D(int argc, char* argv[])
{
  return Test_SpeedToPath_RegularStepGradientDescent_ND<2>(argc, argv);
}

/////////////////////////////////////////////////////////////
// Test for SpeedToPath with RegularStepGradientDescentOptimizer for 3D images
int Test_SpeedToPath_RegularStepGradientDescent_3D(int argc, char* argv[])
{
  return Test_SpeedToPath_RegularStepGradientDescent_ND<3>(argc, argv);
}

/////////////////////////////////////////////////////////////
// Test for SpeedToPath with IterateNeighborhoodOptimizer for 2D images
int Test_SpeedToPath_IterateNeighborhood_2D(int argc, char* argv[])
{
  return Test_SpeedToPath_IterateNeighborhood_ND<2>(argc, argv);
}

/////////////////////////////////////////////////////////////
// Test for SpeedToPath with IterateNeighborhoodOptimizer for 3D images
int Test_SpeedToPath_IterateNeighborhood_3D(int argc, char* argv[])
{
  return Test_SpeedToPath_IterateNeighborhood_ND<3>(argc, argv);
}
