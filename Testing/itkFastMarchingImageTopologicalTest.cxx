#include "itkImage.h"
#include "itkFastMarchingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLabelContourImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"


template <unsigned int VDimension>
int FastMarchingImageFilter( unsigned int argc, char *argv[] )
{
  typedef float InternalPixelType;
  typedef itk::Image<InternalPixelType, VDimension>  InternalImageType;

  typedef unsigned char OutputPixelType;
  typedef itk::Image<OutputPixelType, VDimension> OutputImageType;

  InternalPixelType stoppingValue = atof( argv[5] );

  typedef itk::ImageFileReader< InternalImageType> ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::FastMarchingImageFilter
    <InternalImageType, InternalImageType> FastMarchingFilterType;
  typename FastMarchingFilterType::Pointer fastMarching
    = FastMarchingFilterType::New();
  fastMarching->SetInput( reader->GetOutput() );

  typedef typename FastMarchingFilterType::NodeContainer           NodeContainer;
  typedef typename FastMarchingFilterType::NodeType                NodeType;
  typedef typename FastMarchingFilterType::LabelImageType LabelImageType;

  typedef itk::ImageFileReader<LabelImageType> LabelImageReaderType;
  typename LabelImageReaderType::Pointer labelImageReader = LabelImageReaderType::New();
  labelImageReader->SetFileName( argv[4] );

  try
    {
    labelImageReader->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::LabelContourImageFilter<LabelImageType, LabelImageType> ContourFilterType;
  typename ContourFilterType::Pointer contour = ContourFilterType::New();
  contour->SetInput( labelImageReader->GetOutput() );
  contour->FullyConnectedOff();
  contour->SetBackgroundValue( itk::NumericTraits<typename LabelImageType::PixelType>::Zero );

  try
    {
    contour->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  typename NodeContainer::Pointer alivePoints = NodeContainer::New();
  alivePoints->Initialize();
  unsigned long aliveCount = 0;
  typename NodeContainer::Pointer trialPoints = NodeContainer::New();
  trialPoints->Initialize();
  unsigned long trialCount = 0;

  itk::ImageRegionIteratorWithIndex<LabelImageType> ItL( labelImageReader->GetOutput(),
    labelImageReader->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex<LabelImageType> ItC( contour->GetOutput(),
    contour->GetOutput()->GetLargestPossibleRegion() );
  for( ItL.GoToBegin(), ItC.GoToBegin(); !ItL.IsAtEnd(); ++ItL, ++ItC )
    {
    if( ItC.Get() != itk::NumericTraits<typename LabelImageType::PixelType>::Zero )
      {
      typename LabelImageType::IndexType position = ItC.GetIndex();

      NodeType node;
      const double value = 0.0;

      node.SetValue( value );
      node.SetIndex( position );
      trialPoints->InsertElement( trialCount++, node );
      }
    else if( ItL.Get() != itk::NumericTraits<typename LabelImageType::PixelType>::Zero )
      {
      typename LabelImageType::IndexType position = ItL.GetIndex();

      NodeType node;
      const double value = 0.0;

      node.SetValue( value );
      node.SetIndex( position );
      alivePoints->InsertElement( aliveCount++, node );
      }
    }
  fastMarching->SetTrialPoints(  trialPoints  );
  fastMarching->SetAlivePoints(  alivePoints  );

  fastMarching->SetStoppingValue( stoppingValue );
  fastMarching->SetTopologyCheck( FastMarchingFilterType::None );
  if( argc > 6 && atoi( argv[6] ) == 1 )
    {
    std::cout << "Strict." << std::endl;
    fastMarching->SetTopologyCheck( FastMarchingFilterType::Strict );
    }
  if( argc > 6 && atoi( argv[6] ) == 2 )
    {
    std::cout << "No handles." << std::endl;
    fastMarching->SetTopologyCheck( FastMarchingFilterType::NoHandles );
    }

  try
    {
    fastMarching->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::BinaryThresholdImageFilter
    <InternalImageType, OutputImageType> ThresholdingFilterType;
  typename ThresholdingFilterType::Pointer thresholder
    = ThresholdingFilterType::New();

  thresholder->SetLowerThreshold( 0.0 );
  thresholder->SetUpperThreshold( stoppingValue );
  thresholder->SetOutsideValue( 0 );
  thresholder->SetInsideValue( 1 );
  thresholder->SetInput( fastMarching->GetOutput() );

  try
    {
    thresholder->Update();
    }
   catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }


  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetInput( thresholder->GetOutput() );
  writer->SetFileName( argv[3] );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  if( argc > 7 )
    {
      {
      std::string filename = std::string( argv[7] ) +
        std::string( "LevelSet.nii.gz" );
      typedef itk::ImageFileWriter<InternalImageType> WriterType;
      typename WriterType::Pointer writer = WriterType::New();
      writer->SetInput( fastMarching->GetOutput() );
      writer->SetFileName( filename.c_str() );

      try
        {
        writer->Update();
        }
      catch( itk::ExceptionObject & excep )
        {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }
      }

      {
      std::string filename = std::string( argv[7] ) +
        std::string( "LabelMap.nii.gz" );
      typedef itk::ImageFileWriter< LabelImageType > LabelImageWriterType;
      typename LabelImageWriterType::Pointer mapWriter = LabelImageWriterType::New();
      mapWriter->SetInput( fastMarching->GetLabelImage() );
      mapWriter->SetFileName( filename.c_str() );

      try
        {
        mapWriter->Update();
        }
      catch( itk::ExceptionObject & excep )
        {
        std::cerr << "Exception caught !" << std::endl;
        std::cerr << excep << std::endl;
        return EXIT_FAILURE;
        }
      }
    }

  return EXIT_SUCCESS;
}

int main( int argc, char *argv[] )
{
  if( argc < 6 )
    {
    std::cerr << "Usage: " << argv[0] << " imageDimension";
    std::cerr << " speedImage outputImage seedImage ";
    std::cerr << " stoppingValue [checkTopology] [otherFilePrefix]"
      << std::endl;
    return EXIT_FAILURE;
    }

  switch( atoi( argv[1] ) )
   {
   case 2:
     return FastMarchingImageFilter<2>( argc, argv );
   case 3:
     return FastMarchingImageFilter<3>( argc, argv );
   default:
     std::cerr << "Unsupported dimension" << std::endl;
     return EXIT_FAILURE;
   }
}

