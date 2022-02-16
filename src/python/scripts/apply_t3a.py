"""
Script to apply transform 3a (rigid transform to align images before nonlinear transform based 
registration).

Created by Mukund on 2022-02-14

Usage:

References:
- Reading images - https://simpleitk.readthedocs.io/en/master/IO.html
- Transforming - https://github.com/InsightSoftwareConsortium/SimpleITK-Notebooks/blob/master/Python/21_Transforms_and_Resampling.ipynbs
- Transforming (with output) - http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/21_Transforms_and_Resampling.html
- Transform boundary - http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/Python_html/21_Transforms_and_Resampling.html
- Usefule definitions https://simpleitk.org/SPIE2018_COURSE/images_and_resampling.pdf

"""
import SimpleITK as sitk
import numpy as np


def resample(image, transform, ref_image):
    # Output image Origin, Spacing, Size, Direction are taken from the reference
    # image in this call to Resample
    reference_image = ref_image
    interpolator = sitk.sitkCosineWindowedSinc
    default_value = 100.0
    return sitk.Resample(image, reference_image, transform, interpolator, default_value)


inputImageFileName = "/Users/mraj/Desktop/temp/smaller.png"
# read img to be transformed
reader = sitk.ImageFileReader()
reader.SetImageIO("PNGImageIO")
reader.SetFileName(inputImageFileName)
image = reader.Execute();

print(image)

nda = sitk.GetArrayFromImage(image)
nda = np.moveaxis(nda, [0, 1, 2], [1, 2, 0])
print("shape", nda.shape)
img3d = sitk.GetImageFromArray(nda, isVector=False)

print('img3d', img3d)
# read corresponding transform

tfm = sitk.ReadTransform('/Users/mraj/Desktop/temp/TT.txt')
print(tfm)

# testImage = sitk.Image([10,11,12], sitk.sitkVectorFloat32, 5)
# testImage.SetOrigin((3.0, 14.0, 0))
# testImage.SetSpacing((0.5, 2, 2))

# print(testImage)


# refImageFileName = "/Users/mraj/Desktop/temp/79_01.vsi - 20x.tif"
# reader = sitk.ImageFileReader()
# reader.SetImageIO("TIFFImageIO")
# reader.SetFileName(refImageFileName)
# refImage = reader.Execute();

# nda_ref = sitk.GetArrayFromImage(refImage)
# refImage = sitk.GetImageFromArray(nda_ref, isVector=False)

# extreme_points = [logo.TransformIndexToPhysicalPoint((0,0)), 
#                   logo.TransformIndexToPhysicalPoint((logo.GetWidth(),0)),
#                   logo.TransformIndexToPhysicalPoint((logo.GetWidth(),logo.GetHeight())),
#                   logo.TransformIndexToPhysicalPoint((0,logo.GetHeight()))]


####
translation = sitk.TranslationTransform(3)
translation.SetParameters((-200, -100, 0))
print(translation)
####

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

print("new bounds ", min_x, min_y, max_x, max_y)


print(tfm.TransformPoint(zpt))
# apply transform to images
resampled = resample(img3d, translation.GetInverse(), img3d)
print('width ', img3d.GetWidth())
resampled = resampled_image

# # store image
nda = sitk.GetArrayFromImage(resampled)
nda = np.moveaxis(nda, [0, 1, 2], [2, 0, 1])
print('shape', nda.shape)
img = sitk.GetImageFromArray(nda, isVector=True)
resampled = img

# print ('resampled', resampled)

outputImageFileName = "/Users/mraj/Desktop/temp/smaller_transformed.png"
writer = sitk.ImageFileWriter()
writer.SetFileName(outputImageFileName)
writer.Execute(resampled)










# dimension = 2        
# offset = [2]*dimension # use a Python trick to create the offset list based on the dimension
# translation = sitk.TranslationTransform(dimension, offset)
# print(translation)

# affine = sitk.AffineTransform(2)

# print(affine)



