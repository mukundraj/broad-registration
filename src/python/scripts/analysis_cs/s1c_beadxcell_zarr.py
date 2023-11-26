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
    inp: path to scZarr file created by analysis_sc/s1a_gen_sstab_data_v2.py
    out: path to cell_to_clade_and_class.csv 
    out: path to output puckwise folder (including zarr files)


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
    /v3/s2/wireframes_improved_trans \
    /single_cell/s1/scZarr_321017.zarr \
    /single_cell/s1/cell_to_clade_and_class.csv \
    /cell_spatial/s1/cellspatial_data/cellscores_cshl_231124/ \

Supplementary:

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/cell_spatial/s1/cellspatial_data/cellscores gs://bcdportaldata/cellspatial_data/cellscores

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/cell_spatial/s1/cellspatial_data/cellscores_cshl gs://bcdportaldata/batch_YYMMDD/cellspatial_data/cellscores

gsutil -m rsync -r ~/Desktop/work/data/mouse_atlas/cell_spatial/s1/cellspatial_data/cellscores_cshl_231124 gs://bcdportaldata/batch_231112/cellspatial_data/cellscores_cshl_231124

cd ~/Desktop/work/data/mouse_atlas/cell_spatial/s1/cellspatial_data/cellscores_cshl_231124
fd -p cladeOptions -x  gsutil cp -r {} gs://bcdportaldata/batch_231112/cellspatial_data/cellscores_cshl_231124/{}

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
ip_folder_scZarr = f'{data_root}/{sys.argv[7]}'
op_folder_cell_map = f'{data_root}/{sys.argv[8]}'
op_folder = f'{data_root}/{sys.argv[9]}'

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
    labels_csv_file = f'{label_data_folder}/allen_anno_data_{nis_id_str}.csv'
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

    # write of cladeOptions json file
    cladeOptions_json_file = f'{puck_folder}/cladeOptions.json'
    with open(cladeOptions_json_file, 'w') as outfile:
        json.dump({'cladeOptions':list(unique_clades)}, outfile, separators=(',', ':'))

    # write of cellclassOptions json file
    cellclassOptions_json_file = f'{puck_folder}/cellclassOptions.json'
    with open(cellclassOptions_json_file, 'w') as outfile:
        json.dump({'cellclassOptions':list(unique_cellclasses)}, outfile, separators=(',', ':'))

    # read ccindices json file
    ccindices_json_file = f'{op_folder}/ccindices.json'
    with open(ccindices_json_file) as json_file:
        cc_indices = json.load(json_file)

    # write out zarr with cell score values
    zarr_filename = f'{puck_folder}/cellxbead.zarr'
    store = zarr.DirectoryStore(zarr_filename) # https://zarr.readthedocs.io/en/stable/tutorial.html#storage-alternatives
    nCells, nBeads = counts_X.shape
    dprint('nCells nBeads', nCells, nBeads)
    nAggrs = len(cc_indices) # number of clades and cellclasses aka number of aggretations
    dprint(counts_X)
    root = zarr.group(store=store, overwrite=True)
    zarrX = root.zeros('X', shape=(nCells+nAggrs, nBeads), chunks=(1, nBeads), dtype='f4')
    zarrX[:nCells,:] = np.asarray(counts_X.todense())
    dprint(counts_X.todense())
    maxScoresGroup = root.create_group('maxScores', overwrite=True)
    maxScoresX = maxScoresGroup.zeros('X', shape=(1, nCells+nAggrs), chunks=(1, nCells), dtype='f4')
    # metadata_groupX[:] = np.array([0.0]*nCells);

    # loop over cells and update clade and cellclass contributions from each cell
    tmp_clade_cell_mat = np.zeros((nAggrs, nBeads))
    for cell_idx, cell in enumerate(cells):
        cell_row = counts_X.getrow(cell_idx)
        cell_row_dense = np.squeeze(np.array(cell_row.todense())).astype(float)
        clade, cellclass = cell2cc[cell]
        clade_idx = cc_indices[clade]
        cellclass_idx = cc_indices[cellclass]
        tmp_clade_cell_mat[clade_idx, :] += cell_row_dense
        tmp_clade_cell_mat[cellclass_idx, :] += cell_row_dense

    zarrX[nCells:, :] = tmp_clade_cell_mat

        # zarrX[nCells+clade_idx, :] += cell_row_dense
        # zarrX[nCells+cellclass_idx, :] += cell_row_dense


    # record score sums for each clade and cellclass in current puck
    # cc_score_sums = np.zeros(len(cc_indices)) # write out for each puckid

    cc_score_sums = np.sum(zarrX[nCells:, :], axis=1)

    # make a tmp dir in parent of parent of puck folder if it doesn't exist
    tmp_dir = f'{os.path.dirname(puck_folder)}/../tmp'
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    # write out tmp file with score sums
    tmp_cc_score_sums_file = f'{tmp_dir}/cc_score_sums_{pid}.npy'
    np.save(tmp_cc_score_sums_file, cc_score_sums)

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

    # get maxScores for clades and cellclasses
    for idx, (cc,cc_idx) in enumerate(cc_indices.items()):
        # maxScoresX[0, nCells+cc_idx] = np.max(zarrX[nCells+cc_idx, :])
        maxScoresX[0, nCells+cc_idx] = np.max(tmp_clade_cell_mat[cc_idx, :])


# get cell to clade and class mapping
cell2cc = {}

z = zarr.open(ip_folder_scZarr, mode='r')
cells =  z.obs.clusters[:]
clades = z.metadata.clades[:]
cellclasses = z.metadata.cellclasses[:]

unique_clades = np.unique(clades)
unique_cellclasses = np.unique(cellclasses)


for cell,clade,cellclass in zip(cells,clades,cellclasses):
    cellname = cell.split('=')[1]
    cell2cc[cellname] = [clade, cellclass]

# create a list of clades and classes
uniq_clades_and_classes = unique_clades.tolist() + unique_cellclasses.tolist()

cc_indices = {} # clade and class indices
for idx,cc in enumerate(uniq_clades_and_classes):
    cc_indices[cc] = idx

# write out ccindices json file
ccindices_json_file = f'{op_folder}/ccindices.json'
with open(ccindices_json_file, 'w') as outfile:
    json.dump(cc_indices, outfile, separators=(',', ':'))

pids = list(range(1,208,2))
# pids = list(range(1,5,2))
if __name__ == '__main__':
    start = time.time()
    with Pool(1) as p:
        p.map(process_pid, pids)
    end = time.time()

    # after all pucks are done, process score sums in tmp files to get puckids
    # with max score sum for each clade and cellclass

    # combine tmp files
    tmp_dir = f'{os.path.dirname(op_folder)}/../tmp'

    # get all tmp files
    tmp_files = os.listdir(tmp_dir)
    tmp_files = [f for f in tmp_files if f.startswith('cc_score_sums_')]
    tmp_files = [f for f in tmp_files if f.endswith('.npy')]
    tmp_files = [f'{tmp_dir}/{f}' for f in tmp_files]

    # combine tmp files
    cc_score_sum_mat = np.zeros((len(cc_indices), len(pids)))
    for i,tmp_file in enumerate(tmp_files):
        tmp_arr = np.load(tmp_file)
        cc_score_sum_mat[:, i] = tmp_arr

    # for each clade and cellclass, get puckid with max score sum
    puckinds = np.argmax(cc_score_sum_mat, axis=1)

    # get puckids
    puckids = [pids[i] for i in puckinds]

    # write out puckids and cc_indices in a aggr_info json file
    aggr_info_json_file = f'{op_folder}/aggr_info.json'
    with open(aggr_info_json_file, 'w') as outfile:
        json.dump({'puckids':puckids, 'cc_indices':cc_indices}, outfile, separators=(',', ':'))


    dprint('puckids', puckids)
    dprint(f'Total time {end - start} seconds.')
