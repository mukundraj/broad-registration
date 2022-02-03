// Reads in 3 csv files, 1 NRRD file and outputs 1 csv
// Inputs - From positions, to positions, positions to be transformed and id'ed,
//          and 1 NRRD file with label ids
// Output - a list of label ids stored in a csv 


// References
// - https://itk.org/Doxygen/html/Examples_2RegistrationITKv4_2ThinPlateSplineWarp_8cxx-example.html

#include "itkImage.h"
#include <iostream>
#include "itkTransformFileReader.h"
#include "itkTransformFactoryBase.h"
#include "itkTransformFactory.h"
#include "itkThinPlateSplineKernelTransform.h"

#include "itkPointSet.h"

#include "csv.hpp"


  int
main(int argc, char * argv[])
{
  using namespace csv;

  using ImageType = itk::Image<unsigned short, 3>;
  ImageType::Pointer image = ImageType::New();


  // setup
  constexpr unsigned int ImageDimension = 2;
  using CoordinateRepType = double;
  using TransformType =
    itk::ThinPlateSplineKernelTransform<CoordinateRepType, ImageDimension>;
  using PointType = itk::Point<CoordinateRepType, ImageDimension>;
  using PointSetType = TransformType::PointSetType;
  using PointIdType = PointSetType::PointIdentifier;

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
  CSVReader reader_F("../input/F.csv");  
  CSVReader reader_T("../input/T.csv");  
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

  // Read in query points as csv


  // Compute the TPS transform
  


  // Transform query positions using the transform
  p1[0] = 346.045;
  p1[1] = 455.206;
  p2 = tps->TransformPoint(p1);

  
  std::cout<<"\n";
  std::cout<<p1[0]<<" "<<p1[1]<<" "<<p2[0]<<" "<<p2[1]<<"\n";


  // read NRRD file

  // Record label of each transformed position


  // Write out label of each transformed position

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
