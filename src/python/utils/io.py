"""
IO helper functions

Created by Mukund on 2022-03-22

"""

import csv
import anndata as ann
from scipy.sparse import csr_matrix, hstack
import os
import shutil
import numpy as np
from produtils import dprint

"""
Recreates and empty folder corresponding to given p(uck) id after deleting
existing folder if present.  Returns the path to the created folder.

Created by Mukund on 2022-05-19
"""
def recreate_puck_folder(op_folder, pid):
    puck_folder = f'{op_folder}/puck{pid}'
    if os.path.exists(puck_folder):
        shutil.rmtree(puck_folder)
    os.mkdir(puck_folder)

    return puck_folder

"""

Filters beads that are out of tissue, writes to file remaining bead coords.
Returns indices of on-tissue beads.

Created by Mukund on 2022-05-19
"""

def read_mask_write_beads(pid, ip_folder_labels, ip_folder_chuck_coords, ip_folder_counts, op_folder_gene_csvs):

    # ip_coords_file  = f'{ip_folder_counts}/ad_coords_{str(pid)}.h5ad'
    # coords = ann.read_h5ad(ip_coords_file)

    # coords_X = csr_matrix(coords.X)
    # coords_dense_np = np.array(coords_X.todense())
    # xs = coords_dense_np[:, 0].astype(int).tolist()
    # ys = coords_dense_np[:, 1].astype(int).tolist()
    # zs = coords_dense_np[:, 2].astype(int).tolist()

    nis_id_str = str(pid).zfill(3)
    # get chuck space img coords
    chuck_sp_img_coords_file = f'{ip_folder_chuck_coords}/chuck_sp_img_coords_{nis_id_str}.csv'
    chuck_sp_img_coords = []
    with open(chuck_sp_img_coords_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            chuck_sp_img_coords.append([int(row[0]), int(row[1])])

    # reading csv for label data

    labels_csv_file = f'{ip_folder_labels}/allen_anno_data_{nis_id_str}.csv'
    dprint(labels_csv_file)
    region_names = []
    out_tissue = []
    with open(labels_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            region_names.append(row[4])
            out_tissue.append(row[8])

    dprint(len(region_names), len(out_tissue))

    puck_folder = f'{op_folder_gene_csvs}/puck{pid}'
    # writing coords tsv
    coords_csv_name = f'{puck_folder}/coords.csv'
    # np.savetxt(coords_csv_name, np.array([xs,ys,zs]).T, fmt='%i', header="x,y,z", comments='', delimiter=",")
    in_tissue_inds = []
    with open(coords_csv_name, 'w') as outfile:
        writer = csv.writer(outfile, delimiter=':')
        writer.writerow(['x', 'y', 'rname'])
        for i, status in enumerate(out_tissue):
            if (status=='F'):
                # writer.writerow([xs[i], ys[i], zs[i], region_names[i]])
                writer.writerow([chuck_sp_img_coords[i][0], chuck_sp_img_coords[i][1], region_names[i]])
                in_tissue_inds.append(i)

    return in_tissue_inds

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

""" Reads in new single cell cluster metadata file and returns a dict mapping.
Cluster metadata was provided on 2023-10-16 via Evan slack.

Created by Mukund on 2023-10-17
"""

def get_additional_cluster_metadata(filename):
    additional_metadata = {}
    with open(filename, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            clustername = row[0].split('\t')[1]
            metadata = row[1]
            additional_metadata[clustername] = metadata

    return additional_metadata

