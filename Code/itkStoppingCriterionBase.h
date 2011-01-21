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

#ifndef __itkStoppingCriterionBase_h
#define __itkStoppingCriterionBase_h

#include "itkObject.h"
#include "itkObjectFactoryBase.h"
#include "itkMacro.h"

namespace itk
{
class StoppingCriterionBase : public Object
{
public:
  typedef StoppingCriterionBase Self;
  typedef Object Superclass;
  typedef SmartPointer< Self > Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  itkTypeMacro(StoppingCriterionBase, Object);

  virtual bool IsSatisfied() const = 0;

protected:
  StoppingCriterionBase();
  virtual ~StoppingCriterionBase();

private:
  StoppingCriterionBase( const Self& );
  void operator = ( const Self& );
};
}
#endif
