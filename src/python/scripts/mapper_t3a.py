"""
Applies the initial alignment transform to a set of points from csv and writes out result.

Created by Mukund on 2022-02-16

Usage: python mapper_t3a transform_file.txt input_file.tsv output_file.tsv
Usage example: python src/python/scripts/mapper_t3a.py /Users/mraj/Desktop/temp/TT.txt ./input/points_list.csv ./input/points_list_output.csv
"""
import SimpleITK as sitk
import numpy as np
import sys
import csv

# read transform

# transform_file = "/Users/mraj/Desktop/temp/TT.txt"
transform_file = sys.argv[1]
tfm = sitk.ReadTransform(transform_file)

# read points
input_points_file = sys.argv[2]

transformed_points = []
with open(input_points_file, newline='\n') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        # print(row['first_name'], row['last_name'])
        row = list(map(float, row))

        ## apply transform to point
        point = [row[0], row[1], 0]
        trans_point = tfm.TransformPoint(point)
        transformed_points.append(trans_point)


## write point to csv

output_points_file = sys.argv[3]

with open(output_points_file, 'w', newline='\n') as csvfile:
    writer = csv.writer(csvfile)
    for row in transformed_points:
        writer.writerow(row)

print("done")


