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