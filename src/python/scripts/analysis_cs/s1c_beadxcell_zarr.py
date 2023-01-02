"""
Creates a zarr file for portal from bead x cell anndata files. Also, place
auxaliary files in the same folder - wireframe, nissl, coords data with region
info, cell options json.

Usage:

python s1c_beadxcell_zarr.py
    i/o: data_root
    inp: path to bead_ccf_labels_allbds with bead metadata in csv format
    inp: path to bead_ccf_coords_allbds with bead coords in chuck space in csv format
    inp: input folder with cell x bead anndata files
    inp: path to nissl images folder
    inp: path to atlas wireframe images folder
    out: path to output zarr file


Usage example:

python src/python/scripts/analysis_cs/s1c_beadxcell_zarr.py \
    /Users/mraj/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s7_annotations/bead_ccf_labels_allbds \
    /data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds \
    /cell_spatial/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix_NEW \
    /data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png \
    /data_v3_nissl_post_qc/s7_annotations/allen_labels_imgs/wireframe_trans_bg \
    /cell_spatial/s1/cellspatial_data/cellscores/ \


python src/python/scripts/analysis_cs/s1c_beadxcell_zarr.py \
    /Users/mraj/Desktop/work/data/mouse_atlas \
    /v3/s2/bead_ccf_labels_allbds \
    /data_v3_nissl_post_qc/s3_registered_ss/chuck_img_coords_allbds \
    /cell_spatial/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix_NEW \
    /data_v3_nissl_post_qc/s0_start_formatted_data/transformed_hz_png \
    /v3/s2/wireframes_trans \
    /cell_spatial/s1/cellspatial_data/cellscores_cshl/ \

Supplementary:

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/cell_spatial/s1/cellspatial_data/cellscores gs://bcdportaldata/cellspatial_data/cellscores

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/cell_spatial/s1/cellspatial_data/cellscores_cshl gs://bcdportaldata/cellspatial_data/cellscores

Created by Mukund on 2022-10-24

"""

import sys
import time
from multiprocessing import Pool
from scipy.sparse import csr_matrix, hstack
from produtils import dprint
import anndata as ann
import os
import shutil
import csv
import subprocess
import json
import zarr
import numpy as np

data_root = sys.argv[1]
label_data_folder = data_root+sys.argv[2]
ip_folder_chuck_coords = data_root+sys.argv[3]
ip_folder_cxb = f'{data_root}/{sys.argv[4]}'
ip_folder_nissl = data_root+sys.argv[5]
ip_folder_atlas = data_root+sys.argv[6]
op_folder = f'{data_root}/{sys.argv[7]}'

def process_pid(pid):

    dprint(f'starting pid {pid}..................')
    apid = pid
    if (pid==5 or pid==77 or pid==167):
        apid = pid - 2 ## adjusted pid todo: modify viewer to not require this adjustment

    # ip_coords_file  = f'{in_folder}/ad_coords_{str(apid)}.h5ad'
    ip_counts_file  = f'{ip_folder_cxb}/dd{str(apid)}_CTMapping.h5ad'
    dprint(ip_counts_file)

    counts = ann.read_h5ad(ip_counts_file)
    # coords = ann.read_h5ad(ip_coords_file)

    counts_X = csr_matrix(counts.X).transpose()
    counts_X = csr_matrix(counts.X)
    dprint(counts_X)

    nis_id_str = str(apid).zfill(3)
    # labels_csv_file = f'{label_data_folder}/allen_anno_data_{nis_id_str}.csv'
    labels_csv_file = f'{label_data_folder}/pid_{apid}.csv'
    # dprint(labels_csv_file)
    region_names = []
    with open(labels_csv_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            region_names.append(row[4])

    dprint('region_names', len(region_names))
    # dprint(len(region_names), len(out_tissue))
    puck_folder = f'{op_folder}/puck{pid}'
    dprint('puck_folder', puck_folder)
    if os.path.exists(puck_folder):
        shutil.rmtree(puck_folder)
    os.mkdir(puck_folder)


    # get chuck space img coords
    chuck_sp_img_coords_file = f'{ip_folder_chuck_coords}/chuck_sp_img_coords_{nis_id_str}.csv'
    chuck_sp_img_coords = []
    with open(chuck_sp_img_coords_file, newline='\n') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            chuck_sp_img_coords.append([int(row[0]), int(row[1])])

    dprint('chuck_sp_img_coords', len(chuck_sp_img_coords))
    # copy nissl and atlas images to puck directory
    from_nis_file = f'{ip_folder_nissl}/nis_{nis_id_str}.png'
    to_nis_file = f'{puck_folder}/nis_{nis_id_str}.png'
    from_atlas_file = f'{ip_folder_atlas}/chuck_sp_wireframe_{nis_id_str}.png'
    to_atlas_file = f'{puck_folder}/chuck_sp_wireframe_{nis_id_str}.png'
    # dprint(from_atlas_file)
    # dprint(to_atlas_file)
    subprocess.run(["cp", from_nis_file, to_nis_file])
    subprocess.run(["cp", from_atlas_file, to_atlas_file])

    # writing coords tsv
    coords_csv_name = f'{puck_folder}/coords.csv'
    # dprint(np.shape(coords_dense_np))
    # dprint(coords_csv_name, pid, apid)
    # np.savetxt(coords_csv_name, np.array([xs,ys,zs]).T, fmt='%i', header="x,y,z", comments='', delimiter=",")

    positions = counts.var[['CCF3D_x', 'CCF3D_y', 'CCF3D_z']].values
    position_inds = counts.var.rowName.values.astype(int)
    dprint(positions.shape, position_inds.shape)

    # write out coords csv
    with open(coords_csv_name, 'w') as outfile:
        writer = csv.writer(outfile, delimiter=':')
        # writer.writerow(['x', 'y', 'z', 'rname'])
        writer.writerow(['x', 'y', 'rname'])
        for cidx in position_inds: # iterate over chosen idxs
                writer.writerow([chuck_sp_img_coords[cidx][0], chuck_sp_img_coords[cidx][1], region_names[cidx]])

    # write out cellOptions json file
    cells = list(counts.obs.index)
    cells = [x.split('=')[1] for x in cells]
    cellOptions_json_file = f'{puck_folder}/cellOptions.json'
    gene_options_dict = {'cellOptions':cells}
    with open(cellOptions_json_file, 'w') as outfile:
        json.dump(gene_options_dict, outfile, separators=(',', ':'))

    # write out zarr with cell score values
    zarr_filename = f'{puck_folder}/cellxbead.zarr'
    store = zarr.DirectoryStore(zarr_filename) # https://zarr.readthedocs.io/en/stable/tutorial.html#storage-alternatives
    nCells, nBeads = counts_X.shape
    dprint('nCells nBeads', nCells, nBeads)
    dprint(counts_X)
    root = zarr.group(store=store, overwrite=True)
    zarrX = root.zeros('X', shape=(nCells, nBeads), chunks=(1, nBeads), dtype='f4')
    zarrX[:] = np.asarray(counts_X.todense())
    dprint(counts_X.todense())
    maxScoresGroup = root.create_group('maxScores', overwrite=True)
    maxScoresX = maxScoresGroup.zeros('X', shape=(1, nCells), chunks=(1, nCells), dtype='f4')
    # metadata_groupX[:] = np.array([0.0]*nCells);

    # write out metadata json files
    for cell_idx, cell in enumerate(cells):
        # if (gene_idx%500==0):
        #     collected = gc.collect()
        #     dprint('gene_idx', gene_idx, 'pid', pid, 'collected', collected)
        specific_gene_cnts = counts_X.getrow(cell_idx)
        spec_gene_cnts_dense = np.squeeze(np.array(specific_gene_cnts.todense())).astype(float)
        # spec_gene_cnts_dense = spec_gene_cnts_dense[in_tissue_inds]
        # dprint(np.max(spec_gene_cnts_dense))
        # cell_metadata={"maxCount":np.max(spec_gene_cnts_dense)}
        maxScoresX[0, cell_idx] = np.max(spec_gene_cnts_dense)
        # dprint(cell_idx, maxCountsX[0, cell_idx])
        # gene_cnts=spec_gene_cnts_dense
        # gene_csv_name = f'{puck_folder}/gene_{gene}.csv'
        # np.savetxt(gene_csv_name, gene_cnts, fmt='%i', header="count", comments='',delimiter=',')

        # metadata_json_file = f'{puck_folder}/metadata_cell_{cell}.json'
        # with open(metadata_json_file, 'w') as outfile:
        #     tmp_dict = {'maxCount':str(cell_metadata['maxCount'])}
        #     json.dump(tmp_dict, outfile, separators=(',', ':'))



pids = list(range(1,208,2))
if __name__ == '__main__':
    start = time.time()
    with Pool(1) as p:
        p.map(process_pid, pids)
    end = time.time()
    dprint(f'Total time {end - start} seconds.')
