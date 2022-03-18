"""
For coordinate normalization

Usage:

python
    topleft_x topleft_y
    botright_x botright_y
    topright_x topright_y
    botleft_x botleft_y
    bead_coords.csv
    output_file.csv

Usage example:

python src/python/scripts/v2/s4a_coord_norm.py \
268 274 5497 5502 5497 274 274 5520 \
/Users/mraj/Desktop/sample-hz/coords_143.csv \
/Users/mraj/Desktop/sample-hz/input.csv

Created by Mukund on 2022-03-18

References:

https://stackoverflow.com/questions/57991865/numpy-finding-the-tranformation-matrix-given-2-sets-of-4-points-in-3d-euclidia

"""

import sys
import numpy as np
import csv

topleft_x = sys.argv[1]
topleft_y = sys.argv[2]

botright_x = sys.argv[3]
botright_y = sys.argv[4]

topright_x = sys.argv[5]
topright_y = sys.argv[6]

botleft_x = sys.argv[7]
botleft_y = sys.argv[8]

ip_file = sys.argv[9]
op_file = sys.argv[10]

print(topleft_x, topleft_y, botright_x, botright_y, topright_x, topright_y, botleft_x, botleft_y)
print(ip_file, op_file )

s = np.array([[topleft_x, topleft_y, 1], 
              [topright_x, topright_y, 1],
              [botleft_x, botleft_y, 1],
              [botright_x, botright_y, 1]], dtype=float).T

d = np.array([[0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]], dtype=float).T

M_rec,resid,rank,sing = np.linalg.lstsq(s.T,d.T)
M_rec = M_rec.T


t = np.array([[671, 575, 1], [2566, 2796, 1]],dtype=float)
t = t.T

# read input
input_pts = []
with open(ip_file, newline='\n') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        # print(row['first_name'], row['last_name'])
        row = list(map(float, row))
        point = [row[0], row[1]]
        input_pts.append(point)

print("Number of positions read: ", len(input_pts))

# format input
N = len(input_pts)
pts = np.array(input_pts)
arrays = [pts, np.ones(N).reshape(N,1)]
pts = np.hstack((pts, np.ones(N).reshape(N,1))).T

# transform input
tfmed_pts = M_rec@pts

tfmed_pts = list(tfmed_pts.T)
# print(tfmed_pts)

# write output
with open(op_file, 'w', newline='\n') as csvfile:
    writer = csv.writer(csvfile)
    for row in tfmed_pts:
        writer.writerow([round(row[0], 5), round(row[1], 5)])
