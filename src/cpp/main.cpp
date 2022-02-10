// Reads in 3 csv files, 1 NRRD file and outputs 1 csv
// Inputs - From positions, to positions, positions to be transformed and id'ed,
//          and 1 NRRD file with label ids
// Output - a list of label ids stored in a csv 


// References
// - https://itk.org/Doxygen/html/Examples_2RegistrationITKv4_2ThinPlateSplineWarp_8cxx-example.html
// - https://itk.org/ITKExamples/src/Core/Transform/ApplyAffineTransformFromHomogeneousMatrixAndResample/Documentation.html

#include "itkImage.h"
#include <iostream>
#include "itkTransformFileReader.h"
#include "itkTransformFactoryBase.h"
#include "itkTransformFactory.h"
#include "itkThinPlateSplineKernelTransform.h"
#include "itkImageFileReader.h"
#include "itkResampleImageFilter.h"
#include "itkImageFileWriter.h"

#include "itkPointSet.h"

#include "csv.hpp"
#include <iostream>
#include <fstream>


  int
main(int argc, char * argv[])
{

  if (argc < 4)
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " landmarksFile inputImage ";
    std::cerr << "DeformedImage " << std::endl;
    std::cerr << "deformationField" << std::endl;
    std::cerr << "numargs "<<argc<<std::endl;
    return EXIT_FAILURE;
  }

  using namespace csv;

  using ImageType = itk::Image<unsigned short, 3>;
  ImageType::Pointer image = ImageType::New();


  // setup
  constexpr unsigned int ImageDimension = 2;
  using CoordinateRepType = double;
  using PixelType = unsigned char;
  using InputImageType = itk::Image<PixelType, ImageDimension>;
  using TransformType =
    itk::ThinPlateSplineKernelTransform<CoordinateRepType, ImageDimension>;
  using PointType = itk::Point<CoordinateRepType, ImageDimension>;
  using PointSetType = TransformType::PointSetType;
  using PointIdType = PointSetType::PointIdentifier;
  using FieldVectorType = itk::Vector<float, ImageDimension>;
  using DisplacementFieldType = itk::Image<FieldVectorType, ImageDimension>;
  using ResamplerType =
    itk::ResampleImageFilter<InputImageType, InputImageType>;
  using InterpolatorType =
    itk::LinearInterpolateImageFunction<InputImageType, double>;
  using DeformedImageWriterType = itk::ImageFileWriter<InputImageType>;


  auto      sourceLandMarks = PointSetType::New();
  auto      targetLandMarks = PointSetType::New();
  PointType p1;
  PointType p2;
  PointSetType::PointsContainer::Pointer sourceLandMarkContainer =
    sourceLandMarks->GetPoints();
  PointSetType::PointsContainer::Pointer targetLandMarkContainer =
    targetLandMarks->GetPoints();

  std::string fileName;
  if (argc == 1) // No arguments were provided
  {
    fileName = "../input/Transform.txt";
  }
  else
  {
    fileName = argv[1];
  }

  // Read "From" and "To" csv files
  CSVReader reader_F(argv[1]);  
  CSVReader reader_T(argv[2]);  
  CSVRow row;

  PointIdType id = itk::NumericTraits<PointIdType>::ZeroValue();
  while (reader_F.read_row(row)) {
    std::cout<<row[0].get<int>()<<" "<<row[1].get<double>()<<" "<<row[2].get<double>()<<"\n";
    p1[0] = row[1].get<double>();
    p1[1] = row[2].get<double>();
    sourceLandMarkContainer->InsertElement(id++, p1);
  }

  std::cout<<"\n";
  id = itk::NumericTraits<PointIdType>::ZeroValue();
  while (reader_T.read_row(row)) {
    std::cout<<row[0].get<int>()<<" "<<row[1].get<double>()<<" "<<row[2].get<double>()<<"\n";
    p1[0] = row[1].get<double>();
    p1[1] = row[2].get<double>();
    targetLandMarkContainer->InsertElement(id++, p1);
  }

  auto tps = TransformType::New();
  tps->SetSourceLandmarks(sourceLandMarks);
  tps->SetTargetLandmarks(targetLandMarks);
  tps->ComputeWMatrix();

  std::cout<<"\n";

  // Read in query points as csv

  std::vector<double> points;
  CSVFormat format;
  format.no_header();
  CSVReader reader_P(argv[3], format);
  while (reader_P.read_row(row)) {
    std::cout<<row[0].get<double>()<<" "<<row[1].get<double>()<<"\n";
    points.push_back(row[0].get<double>());
    points.push_back(row[1].get<double>());
  }

  std::cout<<"points "<< points.size()<<std::endl;


  // Compute the TPS transform
  


  // Transform query positions using the transform
  p1[0] = 346.045;
  p1[1] = 455.206;
  p2 = tps->TransformPoint(p1);

  
  std::cout<<"\n";
  // std::cout<<p1[0]<<" "<<p1[1]<<" "<<p2[0]<<" "<<p2[1]<<"\n";

  // write out transformed points

  for (size_t i=0; i< points.size(); i+=2){
    p1[0] = points[i];
    p1[1] = points[i+1];

    p2 = tps->TransformPoint(p1);

    std::cout<<p2[0]<<" "<<p2[1]<<"\n";

    points[i]= p2[0];
    points[i+1] = p2[1];
  }

  // read NRRD file
  using PixelType = unsigned char;
  constexpr unsigned int ImageDimension3D = 3;
  using InputImageType3D = itk::Image<PixelType, ImageDimension3D>;
  using ReaderType = itk::ImageFileReader<InputImageType3D>;

  auto reader = ReaderType::New();
  reader->SetFileName(argv[4]);

  try
  {
    reader->Update();
  }
  catch (const itk::ExceptionObject & excp)
  {
    std::cerr << "Exception thrown " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  using ImageType3D = itk::Image<PixelType, 3>;
  ImageType3D::Pointer image3D = reader->GetOutput(); 

  ImageType3D::RegionType region = image3D->GetLargestPossibleRegion();

  ImageType3D::SizeType size = region.GetSize();

  std::cout << size << std::endl;

  // Record label of each transformed position

  ImageType3D::IndexType indexInside;
  indexInside[0] = 2154;
  indexInside[1] = 2378;
  indexInside[2] = 0;
  ImageType::PixelType value = image3D->GetPixel(indexInside); 
  std::cout << "pixelVal : "<< value<< std::endl;

  // Write out label of each transformed position

  std::stringstream ss; // Can also use ofstream, etc. 
  auto writer = make_csv_writer(ss);
  for (size_t i=0; i<points.size(); i+=2){
    ImageType::PixelType label = image3D->GetPixel(indexInside);
    writer << std::make_tuple(label);
  }

  std::ofstream myfile (argv[5], std::ios::trunc | std::ios::binary);
  if (myfile.is_open()){
    myfile << ss.rdbuf();

  }
  myfile.close();




/////////////////////////////////////////////////////////////////
  // // read "from" image
  // using PixelType = unsigned char;
  // using InputImageType = itk::Image<PixelType, ImageDimension>;
  // using ReaderType = itk::ImageFileReader<InputImageType>;

  // auto reader = ReaderType::New();
  // reader->SetFileName(argv[2]);
 
  // try
  // {
  //   reader->Update();
  // }
  // catch (const itk::ExceptionObject & excp)
  // {
  //   std::cerr << "Exception thrown " << std::endl;
  //   std::cerr << excp << std::endl;
  //   return EXIT_FAILURE;
  // }



  // // tranform "from" image
  // // Set the resampler params
  // InputImageType::ConstPointer inputImage = reader->GetOutput();
  // auto                         resampler = ResamplerType::New();
  // auto                         interpolator = InterpolatorType::New();
  // resampler->SetInterpolator(interpolator);
  // InputImageType::SpacingType   spacing = inputImage->GetSpacing();
  // InputImageType::PointType     origin = inputImage->GetOrigin();
  // InputImageType::DirectionType direction = inputImage->GetDirection();
  // InputImageType::RegionType    region = inputImage->GetBufferedRegion();
  // InputImageType::SizeType      size = region.GetSize();
 
  // // Software Guide : BeginCodeSnippet
  // resampler->SetOutputSpacing(spacing);
  // resampler->SetOutputDirection(direction);
  // resampler->SetOutputOrigin(origin);
  // resampler->SetSize(size);
  // resampler->SetTransform(tps);
  // // Software Guide : EndCodeSnippet
 
  // resampler->SetOutputStartIndex(region.GetIndex());
  // resampler->SetInput(reader->GetOutput());
 
  // // Set and write deformed image
  // auto deformedImageWriter = DeformedImageWriterType::New();
  // deformedImageWriter->SetInput(resampler->GetOutput());
  // deformedImageWriter->SetFileName(argv[3]);
 
  // try
  // {
  //   deformedImageWriter->Update();
  // }
  // catch (const itk::ExceptionObject & excp)
  // {
  //   std::cerr << "Exception thrown " << std::endl;
  //   std::cerr << excp << std::endl;
  //   return EXIT_FAILURE;
  // }
///////////////////////////////////////////////////////////////////


  // using ThinPlateSplineKernelTransformType = itk::ThinPlateSplineKernelTransform<CoordinateRepType, 3>;
  // itk::TransformFactory<itk::ThinPlateSplineKernelTransform<CoordinateRepType, 3>>::RegisterTransform();
  // itk::TransformFactory<ThinPlateSplineKernelTransformType>::RegisterTransform();


  // Register default transforms
  // itk::TransformFactoryBase::RegisterDefaultTransforms();

  // #if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
  //   itk::TransformFileReaderTemplate<CoordinateRepType>::Pointer reader = itk::TransformFileReaderTemplate<CoordinateRepType>::New();
  // #else
  //   itk::TransformFileReader::Pointer writer = itk::TransformFileReader::New();
  // #endif
  //   reader->SetFileName(fileName);
  //   reader->Update();

  //   // Display the transform
  //   std::cout << *(reader->GetTransformList()->begin()) << std::endl;

  std::cout << "ITK Hello World !" << std::endl;
  return EXIT_SUCCESS;
}