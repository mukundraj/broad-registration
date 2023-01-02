"""Generate wireframe in Chuck space by transforming from subsampled original
nissl space.

Usage:

python s2_gen_wireframe.py \
    inp: data root \
    inp: path to file mapping nissl file name to dims \
    inp: geojsons_path \
    inp: path to nissl histolozee .zee file
    out: path to output wireframes

Example:

python src/python/scripts/v3/s2a_gen_wireframe.py \
    ~/Desktop/work/data/mouse_atlas \
    /data_v3_nissl_post_qc/s2_seg_ids/filenames_map.csv \
    /v3/s1/geojson_tfmed \
    /data_v3_nissl_post_qc/s0_start_formatted_data/hz-project.zee \
    /v3/s2/wireframes


Created by Mukund on 2022-12-21
"""



import os
import sys
import json
import geopandas
import matplotlib.pyplot as plt
from shapely.ops import transform
from shapely.ops import unary_union
from geojson import MultiPolygon, Polygon
from produtils import dprint
from pathlib import Path
import sys
path_root = Path(__file__).parents[4]
sys.path.append(str(path_root))
import src.python.utils.io as io
import src.python.utils.transforms as tfms
import numpy as np


data_root = sys.argv[1]
mapperfile_csv = f'{data_root}{sys.argv[2]}'
geojsons_path = f'{data_root}{sys.argv[3]}'
hz_project_file = f'{data_root}{sys.argv[4]}'
wireframes_path = f'{data_root}{sys.argv[5]}'




img_dims = io.get_img_dimensions(mapperfile_csv)
# iterate over geojsons

def tfm_func(x, y, z=None):

    # start from orig nissl img  coords
    x = x
    y = y

    # to orig nissl frac coords
    x = x / width
    y = y / height


    # to hz world coords
    nissl_bb_left = 0
    nissl_bb_top = 0
    x = -x*0.000507*width
    y = -y*0.000507*height


    # to hz world in chuck space
    nissl_id = pid
    mapper_to_new_filename, mapper_to_old_filename = io.get_filenames_map(mapperfile_csv)
    old_nis_filename = mapper_to_old_filename[nissl_id]
    old_nis_id = old_nis_filename.split("_")
    if(len(old_nis_id)==3):
        old_nis_id = int(old_nis_id[1])
    elif (len(old_nis_id)==2):
        old_nis_id = int(old_nis_id[1].split(".")[0])
    elif nissl_id==201:
        old_nis_id=202
    else:
        assert(False)

    tfm_aff = tfms.get_histolozee_affine_tfm_contructed2(hz_project_ss_file=hz_project_file, nis_idx=str(old_nis_id))

    pt_homo = np.array([x, y, 1])
    pt_tfmed = tfm_aff@pt_homo
    x = pt_tfmed[0]
    y = pt_tfmed[1]


    # to frac coords in chuck space
    width_hz = 2.93578
    height_hz = 2.58504
    top = 0.964691 # TODO figure out reason for sign inversion in histolozee export popup!!
    left = 0.670088

    x = (left-x)/width_hz
    y = (top - y)/height_hz


    # to img coords in chuck space
    x = x*4096
    y = 3606-y*3606

    return x, y

for file in os.listdir(geojsons_path)[:]:
    if file.endswith('.geojson'):
        pid = int(file.split('.')[0].split('_')[-1])
        file = f'{geojsons_path}/{file}'
        print(file)
        data = geopandas.read_file(file)
        # print(data.head())
        # data.head()
        # data.plot(aspect=1) # https://stackoverflow.com/questions/67693007/valueerror-box-aspect-and-fig-aspect-must-be-positive
        # plt.show()
        height = img_dims[pid][1]
        width = img_dims[pid][0]

        acronym = []
        ids = []
        name = []
        geometry = []
        for index, row in data.iterrows():

            acro_prop = row['acronym']
            id_prop = row['id']
            name_prop = row['name']
            geom_prop = row['geometry']

            acronym.append(acro_prop)
            ids.append(id_prop)
            name.append(name_prop)

            if geom_prop == None:
                    geometry.append(MultiPolygon())
            else:
                if geom_prop.type == 'Polygon':
                    # g = transform(lambda x, y, z=None: (x, y), geom_prop)
                    g = transform(tfm_func, geom_prop)
                    geometry.append(g)

                elif geom_prop.type == 'MultiPolygon':
                    polygons = []
                    for g in geom_prop.geoms:
                        g = transform(tfm_func, g)


                        polygons.append(g)

                    # external = [p.exterior.coords[:] for p in polygons]
                    # polygons = [p for p in polygons]
                    multiPolygon = unary_union(polygons)
                    # internal = [p.interiors for p in polygons]
                    # multiPolygon = MultiPolygon(external)
                    # multiPolygon = external
                    # for p in polygons:
                    #     dprint(index, len(p.exterior.coords[:]))

                    geometry.append(multiPolygon)

            # print(index)
            # if (index==2):
            #     print('len', len(geom_prop.geoms))

        d = {'acronym': acronym, 'id': id, 'name': name, 'geometry': geometry}
        gdf = geopandas.GeoDataFrame(d)
        # out_file = f'{geojsons_tfmed_path}/pid_{pid}.geojson'
        # gdf.to_file(out_file, driver='GeoJSON')

        # # save a plot
        # data = geopandas.read_file(out_file)
        # print(data.head())
        # data.head()
        w=56.89 # from images of beads transformed to Chuck space
        h=50.08
        # fig = plt.figure(frameon=False)
        fig = plt.figure()
        fig.set_size_inches(w,h)
        # ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        # ax.set_axis_on()
        fig.add_axes(ax)
        # gdf.plot(aspect=1) # https://stackoverflow.com/questions/67693007/valueerror-box-aspect-and-fig-aspect-must-be-positive
        gdf.plot(aspect='auto', ax=ax)  
        pid_str = str(pid).zfill(3)
        plt.xlim([0, 4096-1])
        plt.ylim([0, 3606-1])
        plt.savefig(f'{wireframes_path}/chuck_sp_wireframe_{pid_str}.png', dpi=72)






