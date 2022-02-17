"""
Script to apply transform 3a (rigid transform to align images before nonlinear transform based 
registration).

Created by Mukund on 2022-02-14

Usage: python apply_t3a.py TRANSFORM_FILE IMAGE_TO_BE_TRANSFORMED TRANSFORMED_OUTPUT_FILE_NAME

Usage example: 

python src/python/scripts/apply_t3a.py \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/transforms/79_t3a.txt \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/slide_seq_imgs_rescaled/ss_79.png \
/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/slide_seq_prealigned/pa_ss_79.png

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
    # print('shape', nda.shape)

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

    atfm = sitk.AffineTransform(tfm)
    M = atfm.GetMatrix()

    lastparam = M[-1]

    # Identity cosine matrix .   
    if (lastparam > 0 ):
        output_direction = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    else:
        output_direction = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0]

    # Minimal x,y coordinates are the new origin.
    output_origin = [min_x, min_y, 0]
    # interpolator = sitk.sitkCosineWindowedSinc
    interpolator = sitk.sitkLinear
    # interpolator = sitk.sitkNearestNeighbor
    # Compute grid size based on the physical size and spacing.
    output_size = [int((max_x-min_x)/output_spacing[0]), int((max_y-min_y)/output_spacing[1]), 4]
    resampled_image = sitk.Resample(img3d, output_size, tfm, interpolator, output_origin, output_spacing, output_direction)

    print("new bounds ", min_x, min_y, max_x, max_y)

    # # store image
    nda = sitk.GetArrayFromImage(resampled_image)
    nda = np.moveaxis(nda, [0, 1, 2], [2, 0, 1])
    # print('shape', nda.shape)
    img = sitk.GetImageFromArray(nda, isVector=True)
    resampled = img
    return resampled

# read corresponding transform
transform_file = "/Users/mraj/Desktop/temp/TT.txt"
transform_file = sys.argv[1]
tfm = sitk.ReadTransform(transform_file)
atfm = sitk.AffineTransform(tfm)


# read img to be transformed
inputImageFileName = "/Users/mraj/Desktop/temp/smaller.png"
inputImageFileName = sys.argv[2]
reader = sitk.ImageFileReader()
reader.SetImageIO("PNGImageIO")
reader.SetFileName(inputImageFileName)
image = reader.Execute();

# resample
resampled = transform_img(tfm, image)

# print ('resampled', resampled)

outputImageFileName = "/Users/mraj/Desktop/temp/smaller_transformed.png"
outputImageFileName = sys.argv[3]
writer = sitk.ImageFileWriter()
writer.SetFileName(outputImageFileName)
writer.Execute(resampled)


