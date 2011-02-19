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

#ifndef __itkArrivalFunctionToPathFilter_txx
#define __itkArrivalFunctionToPathFilter_txx

#include "itkArrivalFunctionToPathFilter.h"

namespace itk
{
template <class TInputImage, class TOutputPath>
ArrivalFunctionToPathFilter<TInputImage,TOutputPath>
::ArrivalFunctionToPathFilter()
{
  m_TerminationValue = 1.0;
  m_CurrentOutput = 0;
}

template <class TInputImage, class TOutputPath>
ArrivalFunctionToPathFilter<TInputImage,TOutputPath>
::~ArrivalFunctionToPathFilter()
{
}

template <class TInputImage, class TOutputPath>
void
ArrivalFunctionToPathFilter<TInputImage,TOutputPath>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if ( this->GetInput() )
    {
    InputImagePointer image =
      const_cast< InputImageType * >( this->GetInput() );
    image->SetRequestedRegionToLargestPossibleRegion();
    }
}

template<class TInputImage, class TOutputPath>
unsigned int
ArrivalFunctionToPathFilter<TInputImage,TOutputPath>
::GetNumberOfPathsToExtract() const
{
  return m_PointList.size();
}

template<class TInputImage, class TOutputPath>
const typename ArrivalFunctionToPathFilter<TInputImage,TOutputPath>::PointType &
ArrivalFunctionToPathFilter<TInputImage,TOutputPath>
::GetNextEndPoint()
{
  return m_PointList[m_CurrentOutput];
}

template<class TInputImage, class TOutputPath>
typename ArrivalFunctionToPathFilter<TInputImage,TOutputPath>::InputImageType *
ArrivalFunctionToPathFilter<TInputImage,TOutputPath>
::ComputeArrivalFunction()
{
  InputImageType * function = (InputImageType *)this->ProcessObject::GetInput(0);
  return function;
}

template <class TInputImage, class TOutputPath>
void
ArrivalFunctionToPathFilter<TInputImage,TOutputPath>
::GenerateData( void )
{
  // Get the input
  InputImagePointer input =
    const_cast< InputImageType * >( this->GetInput() );
  if ( input.IsNull() )
    {
    itkExceptionMacro( "Input image must be provided" );
    return;
    }

  // Check the number of paths is not none
  unsigned int numberOfOutputs = GetNumberOfPathsToExtract();
  if ( numberOfOutputs == 0 )
    {
    itkExceptionMacro( "At least one path must be specified for extraction" );
    }
  this->ProcessObject::SetNumberOfRequiredOutputs( numberOfOutputs );

  // Setup cost function
  if ( m_CostFunction.IsNull() )
    {
    m_CostFunction = CostFunctionType::New();
    }

  // Setup optimizer
  if ( m_Optimizer.IsNull() )
    {
    typename InputImageType::SpacingType spacing = input->GetSpacing();

    // Compute the minimum spacing
    double minspacing = static_cast< double >( spacing[0] );
    double s;

    for (unsigned int dim=1; dim < InputImageDimension; dim++)
      {
      s = static_cast< double >( spacing[dim] );
      if ( s < minspacing )
        {
        minspacing = s;
        }
      }
    // Create default optimizer
    typename DefaultOptimizerType::Pointer defaultOptimizer =
        DefaultOptimizerType::New();
    defaultOptimizer->SetNumberOfIterations( 1000 );
    defaultOptimizer->SetMaximumStepLength( 1.5*minspacing );
    defaultOptimizer->SetMinimumStepLength( 0.5*minspacing );
    defaultOptimizer->SetRelaxationFactor( 0.999 );
    m_Optimizer = defaultOptimizer;
    }
  m_Optimizer->SetCostFunction( m_CostFunction );

  // Observe optimizer
  typename CommandType::Pointer observer = CommandType::New();
  observer->SetFilter( this );
  unsigned long observerTag =
      m_Optimizer->AddObserver( itk::IterationEvent(), observer );

  // Do for each output
  for ( unsigned int n=0; n < numberOfOutputs; n++ )
    {
    // Set the output index
    // NOTE: m_CurrentOutput is used in Execute() and GetNextEndPoint()
    m_CurrentOutput = n;

    // Make the output
    OutputPathPointer output
       = static_cast<TOutputPath*>( this->MakeOutput(n).GetPointer() );
    this->ProcessObject::SetNthOutput( n, output.GetPointer() );

    // Compute the arrival function
    InputImagePointer function = this->ComputeArrivalFunction( );

    if ( m_CostFunction->GetImage() != function )
      {
      m_CostFunction->SetImage( function );
      m_CostFunction->Initialize( );
      }

    // Get the end point to back propagate from
    PointType pointEnd = this->GetNextEndPoint();

    // Convert end point to parameters type
    typename CostFunctionType::ParametersType end(InputImageDimension);
    for (unsigned int i=0; i<InputImageDimension; i++)
      {
      end[i] = pointEnd[i];
      }

    // Initialize optimizer
    m_Optimizer->SetInitialPosition( end );

    // Use optimizer to back propagate from end point
    m_Optimizer->StartOptimization( );
    }

  // Clean up by removing observer
  m_Optimizer->RemoveObserver( observerTag );
}

template <class TInputImage, class TOutputPath>
void
ArrivalFunctionToPathFilter<TInputImage,TOutputPath>
::Execute( const Object * object, const EventObject & )
{
  // Cast object to optmizer
  typename OptimizerType::Pointer optimizer = (OptimizerType*)
      dynamic_cast< const OptimizerType* >( object );
  if ( optimizer.IsNull() )
    {
    return;
    }

  // Get current position and value
  typename OptimizerType::ParametersType
      currentParameters = optimizer->GetCurrentPosition();

  unsigned int lenParameters = currentParameters.GetSize();

  if ( lenParameters != InputImageDimension )
    {
    return;
    }

  typename OptimizerType::MeasureType
      currentValue = optimizer->GetValue( currentParameters );

  // Check the current value is above given termination threshold
  if ( currentValue < m_TerminationValue )
    {
    return;
    }

  // Convert parameters to point
  bool valid = false;
  unsigned int numparams = optimizer->GetCurrentPosition().GetSize();

  PointType point;
  point.Fill( 0.0 );

  for ( unsigned int i=0; i<numparams; i++ )
    {
    point[i] = optimizer->GetCurrentPosition()[i];
    valid = true;
    }
  if ( !valid )
    {
    return;
    }

  // Convert point to continuous index
  InputImagePointer input = const_cast<InputImageType*>( this->GetInput() );
  ContinuousIndexType cindex;
  input->TransformPhysicalPointToContinuousIndex( point, cindex );

  // Add point as vertex in path
  typename OutputPathType::Pointer output = this->GetOutput( m_CurrentOutput );
  output->AddVertex( cindex );
}

template<class TInputImage, class TOutputPath>
void
ArrivalFunctionToPathFilter<TInputImage,TOutputPath>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "TerminationValue: " << m_TerminationValue << std::endl;
  os << indent << "NumberOfEndPoints: " << m_PointList.size() << std::endl;
}


} // end namespace itk

#endif
