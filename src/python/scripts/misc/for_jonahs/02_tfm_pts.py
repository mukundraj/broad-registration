"""
Transform points to chuck space using Slicer's transform file in csv format

python \
    ip: path to coords to transform
    ip: path to transform in csv format
    op: op folder

python src/python/scripts/misc/for_jonahs/02_tfm_pts.py \
     "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/ip_dlgn/dlgn_markers_bead_coords_df.csv" \
     "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/v2/tfm.txt" \
     "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/v2"

Created by Mukund on 2024-09-20
""" 

import numpy as np
import csv
import subprocess
import SimpleITK as sitk
import sys

ip_dlgn = sys.argv[1]
tfm_path = sys.argv[2]

# ip_all = "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/v2/091824_thalamus_all_bead_coords.csv"
# ip_f = "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/v2/From.fcsv"
# ip_t = "/Users/mraj/Desktop/work/data/mouse_atlas/misc/for_jonahs/v2/To.fcsv"
op_folder = sys.argv[3]

# read csv file
data_all = np.genfromtxt(ip_all, delimiter=',', dtype=None, encoding=None)
data_dlgn = np.genfromtxt(ip_dlgn, delimiter=',', dtype=None, encoding=None)

coords = data_dlgn[1:, 1:3]

# convert to a numeric
coords = coords.astype(np.float16)

print(data_dlgn)

min_x = np.min(coords[:, 0])
max_x = np.max(coords[:, 0])

min_y = np.min(coords[:, 1])
max_y = np.max(coords[:, 1])

# make coords positive
coords[:, 0] = coords[:, 0] - min_x
coords[:, 1] = coords[:, 1] - min_y

# get new extents
max_x = np.max(coords[:, 0])
max_y = np.max(coords[:, 1])
min_x = np.min(coords[:, 0])
min_y = np.min(coords[:, 1])

print (max_x, min_x, max_y, min_y)

# read the tfms

def tfm_pts(tfm_path, pts, transform_ycoords): # bool tfm_ycoords

    tfm = sitk.ReadTransform(tfm_path)
    # print(tfm)
    M = str(tfm).split('\n')[10:13]
    M = [i.strip().split(' ') for i in M]
    M = [[float(j) for j in i] for i in M]
    M = np.array(M)[:-1,:-1] # discarding 3rd dimension

    O = str(tfm).split('\n')[13][12:-1]
    O = [float(i) for i in O.split(',')]
    O = np.array(O[:-1]).reshape(2,1)

    C = str(tfm).split('\n')[14][12:-1]
    C = [float(i) for i in C.split(',')]
    C = np.array(C[:-1]).reshape(2,1)

    T = str(tfm).split('\n')[15][17:-1]
    T = [float(i) for i in T.split(',')]
    T = np.array(T[:-1]).reshape(2,1)

    A = np.zeros((3,3))
    A[:2,:2] = M
    A[:2,2] = T[:,0]
    A[2,2] = 1

    pts_tmp = pts[:2, :]

    y_range = np.max(pts_tmp[1, :])

    # subtract y axis from 525 to get correct y axis
    if transform_ycoords == -1:
        pts_tmp[1, :] =   pts_tmp[1, :]
    else:
        pts_tmp[1, :] =  y_range - pts_tmp[1, :]

    pts_tfmed = []
    pts_tmp = np.vstack((pts_tmp, np.zeros((1, pts_tmp.shape[1]))))
    for i in range(pts_tmp.shape[1]):
        pts_tfmed.append(tfm.GetInverse().TransformPoint(pts_tmp[:, i].T))
    pts_tfmed = np.array(pts_tfmed).T


    # # print('M\n', M)
    # # print('O\n', O)
    # # print('C\n', C)
    # # print('T\n', T)
    # # print('A\n', A)

    # # print A
    # # print('A\n', A)
    # # invert A
    # A_inv = np.linalg.inv(A)
    # A_inv = np.array([[-2.09403934e-01, -8.25087080e-01,  5.27251472e+02],
    #                  [ 8.25087080e-01, -2.09403934e-01,  8.05165971e+02],
    #                  [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    # # A_inv = np.array([[-2.09403934e-01, -8.25087080e-01, 0],
    # #                  [ 8.25087080e-01, -2.09403934e-01, 0],
    # #                  [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])

    # # print('A_inv\n', A_inv)

    # pts_tfmed = np.matmul(A_inv, pts)
    return A, pts_tfmed

op_file_tfmed_pts = f'{op_folder}/tfmed.csv'

# coords = np.vstack((np.array([1869.0435965349486,604.2925401958946]), coords))
# coords = np.vstack((np.array([2217.458539526952,799.7448252889699]), coords))

 
A, pts_tfmed = tfm_pts(tfm_path, coords.T, -1)

print(np.shape(coords))
print(np.shape(pts_tfmed))

np.savetxt(op_file_tfmed_pts, pts_tfmed[:2].T, delimiter=",")
# np.savetxt(op_file_tfmed_pts, pts_tfmed, delimiter=",")

