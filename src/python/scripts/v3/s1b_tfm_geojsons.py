"""Transform Partha geojson contours to original nissl space.

Usage:

Example:

python s1b_tfm_geojsons.py \
    inp: data root
    inp: geojsons_path
    inp: id to tiff file path
    inp: tiff dims file path
    inp: acro key map path (from s1_acro_key_map.py)
    inp: orig dims map path (from s1_gen_origdims_map.py)
    out: output path for transformed geojsons
    out: output path for transformed geojsons images

python src/python/scripts/v3/s1b_tfm_geojsons.py \
    ~/Desktop/work/data/mouse_atlas \
    /v3/s0/geojson \
    /data_v3_nissl_post_qc/s2_seg_ids/id_to_tiff_mapper.csv \
    /data_v3_nissl_post_qc/s2_seg_ids/filenames_map.csv \
    /v3/s1/acro_key_map.json \
    /v3/s1/orig_dims_map.csv \
    /v3/s1/geojson_tfmed \
    /v3/s1/geojson_tfmed_imgs \


Created on 2022-12-19
"""


import geopandas
import matplotlib.pyplot as plt
import sys
import os
from produtils import dprint
import shapely
from geojson import MultiPolygon, Polygon
from shapely import wkt
from shapely.geometry import Point
import itertools
from shapely.ops import unary_union
import csv
from shapely.ops import transform
import json

data_root = sys.argv[1]
geojsons_path = f'{data_root}{sys.argv[2]}'
id_to_tiff_file_path = f'{data_root}{sys.argv[3]}'
tiff_dims_file_path = f'{data_root}{sys.argv[4]}'
acro_key_map_path = f'{data_root}{sys.argv[5]}'
orig_dims_map_path = f'{data_root}{sys.argv[6]}'
geojsons_tfmed_path = f'{data_root}{sys.argv[7]}'
geojsons_tfmed_imgs = f'{data_root}{sys.argv[8]}'

id_to_tiff = {}
with open(id_to_tiff_file_path) as f:
    reader = csv.reader(f)
    for row in reader:
        id_to_tiff[row[0]] = row[1]

tiff_dims = {}
with open(tiff_dims_file_path) as f:
    reader = csv.reader(f)
    for row in reader:
        tiff_dims[row[0]] = [int(row[2]), int(row[3])]

# read acro key map
with open(acro_key_map_path) as f:
    acro_key_map = json.load(f)

# read orig dims map
orig_dims_map = {}
with open(orig_dims_map_path) as f:
    reader = csv.reader(f)
    for row in reader:
        orig_dims_map[row[0]] = [int(row[2]), int(row[3])]

# files = os.listdir(geojsons_path)
# files = [f for f in files if f.endswith('.geojson')]
files = []
idxs = []
start_pid = 1
end_pid = 207
for i in range (int((start_pid+1)/2), int((end_pid+1)/2+1)):
    i_str = str(i).zfill(4)
    fname = f'atlas_to_MD743&742-N1-2019.02.13-21.15.52_MD743_1_{i_str}.geojson'
    files.append(fname)
    idxs.append(i)

print(idxs)
for idx, file in zip(idxs, files):
    acronym = []
    id = []
    name = []
    geometry = []

    pid = int(idx*2 - 1)

    if (pid==5 or pid==77 or pid==167):
        continue

    tiff_file = id_to_tiff[str(pid)]
    cur_tiff_dims = tiff_dims[f'nis_{tiff_file}'[:-1]]
    dprint(idx, cur_tiff_dims, tiff_file, 'pid', pid)
    file_full_path = f'{geojsons_path}/{file}'

    data = geopandas.read_file(file_full_path)
    # dprint(data.bounds)
    dprint(data.total_bounds)

    # print(data.head())

    for index, row in data.iterrows():

        acro_prop = row['acronym']
        id_prop = acro_key_map[acro_prop]
        name_prop = row['name']
        geom_prop = row['geometry']
        if (len(geom_prop.geoms)>0):
            # dprint(len(mp.geoms))
            polygons = []
            for g in geom_prop.geoms:
                # rotate
                g = transform(lambda x, y, z=None: (-y, x), g)

                # # convert to img coords
                # g = transform(lambda x, y, z=None: (x, cur_tiff_dims[1]-y), g)


                # translate
                full_tiff_dims = [0.3428*orig_dims_map[str(pid)][0], 0.3428*orig_dims_map[str(pid)][1]]

                dx = (full_tiff_dims[0] ) * 0.5 # this needs to be in microns
                dy = (full_tiff_dims[1] ) * 0.5
                # g = transform(lambda x, y, z=None: (cur_tiff_dims[0] - (x+dx), y-dy), g)
                g = transform(lambda x, y, z=None: ((x+dx), y+dy), g)

                # scale to image space (um to pixels) multiply by pixels per um
                res = 0.3428
                g = transform(lambda x, y, z=None: (x/res, y/res), g)

                # scale to subsampled space
                g = transform(lambda x, y, z=None: (x/8, y/8), g)

                polygons.append(g)

            # external = [p.exterior.coords[:] for p in polygons]
            # polygons = [p for p in polygons]
            multiPolygon = unary_union(polygons)
            # internal = [p.interiors for p in polygons]
            # multiPolygon = MultiPolygon(external)
            # multiPolygon = external
            # for p in polygons:
            #     dprint(index, len(p.exterior.coords[:]))

            acronym.append(acro_prop)
            id.append(id_prop)
            name.append(name_prop)
            geometry.append(multiPolygon)
        else:
            acronym.append(acro_prop)
            id.append(id_prop)
            name.append(name_prop)
            geometry.append(MultiPolygon())

        # print(index)
        # if (index==2):
        #     print('len', len(geom_prop.geoms))

    d = {'acronym': acronym, 'id': id, 'name': name, 'geometry': geometry}
    gdf = geopandas.GeoDataFrame(d)
    out_file = f'{geojsons_tfmed_path}/pid_{pid}.geojson'
    gdf.to_file(out_file, driver='GeoJSON')

    # save a plot
    data = geopandas.read_file(out_file)
    # print(data.head())
    # data.head()
    data.plot(aspect=1) # https://stackoverflow.com/questions/67693007/valueerror-box-aspect-and-fig-aspect-must-be-positive

    pid_str = str(pid).zfill(3)
    file = f'/Users/mraj/Desktop/work/data/mouse_atlas/data_v3_nissl_post_qc/s4_bead_to_segid/bead_to_nis_coords/bead_to_nis_coords_{pid_str}.csv'
    # read csv file
    with open(file) as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            x = float(row[1])
            # y = cur_tiff_dims[1] - float(row[2])
            y = float(row[2])
            if i%50==0:
                plt.scatter(x, y, s=1, c='r')
    plt.savefig(f'{geojsons_tfmed_imgs}/pid_{pid}.png')


    dprint ('len geom', len(geometry), 'len acronym', len(acronym))





# d = {'acronym':['3', 'sdf'], 'name':[3,3], 'id': [1, 2], 'geometry': [Point(1, 2), Point(2, 1)]}
# gdf = geopandas.GeoDataFrame(d)
# out_file = f'{geojsons_tfmed_path}/test.geojson'
# gdf.to_file(out_file, driver='GeoJSON')
# print (gdf)


exit(0)

for index, row in data.iterrows():
    # print(row)
    # print(row['geometry'])
    mp = row['geometry']
    if (len(mp.geoms)>0):
        # dprint(len(mp.geoms))
        p = mp.geoms[0]
        p2 = shapely.transform(p, lambda x: x+1, lambda y: y+1)
        print(type(p2))
        # print(p2)
        # p = Polygon(p)
        # mp2 = MultiPolygon([p])
        # mp2 = MultiPolygon([wkt.loads(pw)])
        # ob = MultiPolygon([
        #     (
        #         ((0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)),
        #         [((0.1,0.1), (0.1,0.2), (0.2,0.2), (0.2,0.1))]
        #     )
        # ])
        mp2 = MultiPolygon([p.exterior.coords[:], p2.exterior.coords[:]])
        dprint(mp2)

        
        exit(0)


