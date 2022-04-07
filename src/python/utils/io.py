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