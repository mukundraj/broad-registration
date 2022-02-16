"""
Script to apply transform 3a (rigid transform to align images before nonlinear transform based 
registration).

Created by Mukund on 2022-02-14

Usage: python apply_t3a.py TRANSFORM_FILE IMAGE_TO_BE_TRANSFORMED

References:
- Reading images - https://simpleitk.readthedocs.io/en/master/IO.html
- Transforming - https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/21_Transforms_and_Resampling.ipynbs
- Transforming (with output) - http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/21_Transforms_and_Resampling.html
- Transform boundary - http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/21_Transforms_and_Resampling.html
- Usefule definitions https://simpleitk.org/SPIE2018_COURSE/images_and_resampling.pdf

"""
import SimpleITK as sitk
import numpy as np
import sys


def transform_img(tfm, image):

    # convert to 3D image
    nda = sitk.GetArrayFromImage(image)
    nda = np.moveaxis(nda, [0, 1, 2], [1, 2, 0])
    img3d = sitk.GetImageFromArray(nda, isVector=False)


    zpt = img3d.TransformIndexToPhysicalPoint((0, 0, 0))
    extreme_points = [img3d.TransformIndexToPhysicalPoint((0,0,0)), 
                      img3d.TransformIndexToPhysicalPoint((img3d.GetWidth(),0, 0)),
                      img3d.TransformIndexToPhysicalPoint((img3d.GetWidth(),img3d.GetHeight(), 0)),
                      img3d.TransformIndexToPhysicalPoint((0,img3d.GetHeight(), 0))]
    extreme_points_transformed = [tfm.GetInverse().TransformPoint(pnt) for pnt in extreme_points]
    min_x = min(extreme_points_transformed)[0]
    min_y = min(extreme_points_transformed, key=lambda p: p[1])[1]
    max_x = max(extreme_points_transformed)[0]
    max_y = max(extreme_points_transformed, key=lambda p: p[1])[1]

    # Use the original spacing (arbitrary decision).
    output_spacing = img3d.GetSpacing()
    # Identity cosine matrix (arbitrary decision).   
    output_direction = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    # Minimal x,y coordinates are the new origin.
    output_origin = [min_x, min_y, 0]
    # Compute grid size based on the physical size and spacing.
    output_size = [int((max_x-min_x)/output_spacing[0]), int((max_y-min_y)/output_spacing[1]), 4]
    resampled_image = sitk.Resample(img3d, output_size, tfm, sitk.sitkLinear, output_origin, output_spacing, output_direction)

    # print("new bounds ", min_x, min_y, max_x, max_y)

    # # store image
    nda = sitk.GetArrayFromImage(resampled_image)
    nda = np.moveaxis(nda, [0, 1, 2], [2, 0, 1])
    # print('shape', nda.shape)
    img = sitk.GetImageFromArray(nda, isVector=True)
    resampled = img
    return resampled

inputImageFileName = "/Users/mraj/Desktop/temp/smaller.png"
inputImageFileName = sys.argv[2]
# read img to be transformed
reader = sitk.ImageFileReader()
reader.SetImageIO("PNGImageIO")
reader.SetFileName(inputImageFileName)
image = reader.Execute();

# read corresponding transform

tranform_file = "/Users/mraj/Desktop/temp/TT.txt"
transform_file = sys.argv[1]
tfm = sitk.ReadTransform(transform_file)

resampled = transform_img(tfm, image)

# print ('resampled', resampled)

outputImageFileName = "/Users/mraj/Desktop/temp/smaller_transformed.png"
writer = sitk.ImageFileWriter()
writer.SetFileName(outputImageFileName)
writer.Execute(resampled)


