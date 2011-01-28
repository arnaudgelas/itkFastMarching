/*=========================================================================
  
  Project:   MinimalPath
  Program:   Insight Segmentation & Registration Toolkit
  Module:    MinimalPathTestDriver.cxx
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

//This file is the test driver for the ComposeRGBA project.
//It is responsible for registering all test methods.

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#include <iostream>
#include "itkTestMain.h"

void RegisterTests()
{
	//Register test methods  
	REGISTER_TEST(Test_SpeedToPath_GradientDescent_2D);
	REGISTER_TEST(Test_SpeedToPath_GradientDescent_3D);
	REGISTER_TEST(Test_SpeedToPath_RegularStepGradientDescent_2D);
	REGISTER_TEST(Test_SpeedToPath_RegularStepGradientDescent_3D);
	REGISTER_TEST(Test_SpeedToPath_IterateNeighborhood_2D);
	REGISTER_TEST(Test_SpeedToPath_IterateNeighborhood_3D);
}
