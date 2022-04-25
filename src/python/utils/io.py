"""
IO helper functions

Created by Mukund on 2022-03-22

"""

import csv

"""
Reads in path to corners csv and returns a map from nissl_id to corners info.

Created by Mukund on 2022-03-22
"""
def get_corners_info(corners_csv):

    corners = {}
    with open(corners_csv, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            # print(row['first_name'], row['last_name'])
            row = list(map(int, row))
            corners[row[0]] = {'topleft_x':row[1], 'topleft_y':row[2],
                               'botright_x':row[3], 'botright_y':row[4],
                               'topright_x':row[5], 'topright_y':row[6],
                               'botleft_x':row[7], 'botleft_y':row[8]}


    return corners

""" Reads in filename mapper file to convert from between old to new nissl
filenames(and vice versa) and returns two dicts.

Created by Mukund on 2022-04-20
"""
def get_filenames_map(mapperfile_csv):
    map_to_new_filename = {}
    map_to_old_filename = {}
    with open(mapperfile_csv, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        # next(reader)
        for row in reader:
            nissl_id = int(row[1])
            if (nissl_id>0):
                nis_id_str = str(nissl_id).zfill(3)
                map_to_new_filename[row[0]] = nis_id_str
                map_to_old_filename[nissl_id] = row[0]

    return map_to_new_filename, map_to_old_filename


""" Reads in filename mapper file and returns dict containing img dimensions.
The mapper file is same file used to map between old and new nissl names, with
additional columns added later storing corresponding image dimensions.

Created by Mukund on 2022-04-20
"""
def get_img_dimensions(mapperfile_csv):

    img_dims = {}
    with open(mapperfile_csv, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        # next(reader)
        for row in reader:
            nissl_id = int(row[1])
            if (nissl_id>0):
                img_dims[nissl_id] = [int(row[2]), int(row[3])]

    return img_dims
