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

#ifndef __itkIterateNeighborhoodOptimizer_h
#define __itkIterateNeighborhoodOptimizer_h

#include "itkArray.h"
#include "itkSingleValuedNonLinearOptimizer.h"

namespace itk
{

/** \class IterateNeighborhoodOptimizer
 * \brief Finds the local minima/maxima by iteratively choosing the
 *        minimum/maximum value in a neighborhood.
 *
 * This optimizer is designed to operate on a monotonic cost function
 * WITHOUT using gradient information (derivatives). The user must set
 * the Neighborhood size, and optionally the connectivity.
 *
 * \ingroup Numerics Optimizers
 */
class ITK_EXPORT IterateNeighborhoodOptimizer :
    public SingleValuedNonLinearOptimizer
{
public:
  /** Standard class typedefs. */
  typedef IterateNeighborhoodOptimizer   Self;
  typedef SingleValuedNonLinearOptimizer Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( IterateNeighborhoodOptimizer, SingleValuedNonLinearOptimizer );

  /** Configure whether the local maxima or minima is found.
   *  The default is to minimize the cost function (maximize = false ).*/
  itkGetConstReferenceMacro( Maximize, bool );
  itkSetMacro( Maximize, bool );
  itkBooleanMacro( Maximize );
  bool GetMinimize( ) const
    { return !m_Maximize; }
  void SetMinimize(bool v)
    { this->SetMaximize(!v); }
  void MinimizeOn()
    { this->MaximizeOff(); }
  void MinimizeOff()
    { this->MaximizeOn(); }

  /** Advance one step. */
  virtual void AdvanceOneStep( void );

  /** Start optimization. */
  void StartOptimization( void );

  /** Resume previously stopped optimization with current parameters
   * \sa StopOptimization. */
  void ResumeOptimization( void );

  /** Stop optimization.
   * \sa ResumeOptimization */
  void StopOptimization( void );

  /**
   * Get/set whether the nieghborhood is defined by face connectivity or
   * by face+edge+vertex connectivity.  Default is FullyConnectedOn.
   */
  itkSetMacro(FullyConnected, bool);
  itkGetConstReferenceMacro(FullyConnected, bool);
  itkBooleanMacro(FullyConnected);

  /** Get/set the nieghborhood size (in physical space).
   *  The default is [1.0,1.0] and MUST be specified for all 3-D images
   *  and 2-D images with non-unity spacing. */
  typedef Array<double> NeighborhoodSizeType;
  itkSetMacro( NeighborhoodSize, NeighborhoodSizeType );
  itkGetConstReferenceMacro( NeighborhoodSize, NeighborhoodSizeType );

  /** Get the current iteration number. */
  itkGetConstMacro( CurrentIteration, unsigned int );

  /** Get the current value. */
  itkGetConstReferenceMacro( CurrentValue, double );

protected:
  IterateNeighborhoodOptimizer();
  virtual ~IterateNeighborhoodOptimizer() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  IterateNeighborhoodOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  bool                 m_Stop;
  bool                 m_Maximize;
  bool                 m_FullyConnected;
  double               m_CurrentValue;
  unsigned long        m_CurrentIteration;
  NeighborhoodSizeType m_NeighborhoodSize;

};

} // end namespace itk

#endif
