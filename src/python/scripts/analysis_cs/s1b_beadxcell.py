"""
Read and process bead x cell score matrices. Generate csv with positions for
paraview visualization.


Usage:

python s1b_beadxcell.py
    out: output folder for csv with positions

Usage example:

python src/python/scripts/analysis_sc/s1b_beadxcell.py \
    /Users/mraj/Desktop/work/data/mouse_atlas/single_cell/s1/for3d \


Supplementary:

gsutil -m rsync -r gs://macosko_data/jlanglie/scp/03_All_MBASS_Mapping_Mega_Matrix ~/Desktop/work/data/mouse_atlas/single_cell/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix

Created by Mukund on 2022-10-13

"""

from matplotlib.cbook import file_requires_unicode
from produtils import dprint
import anndata as ann
import numpy as np
import sys
import vtk

in_folder = "/Users/mraj/Desktop/work/data/mouse_atlas/single_cell/s0/raw_beadxctype/03_All_MBASS_Mapping_Mega_Matrix"

inp_file = f'{in_folder}/dd1_CTMapping.h5ad'
out_folder = sys.argv[1]

out_file = f'{out_folder}/fvals.csv'

data = ann.read_h5ad(inp_file)

# dprint(data)
# dprint(data.X)
# dprint(data.obs)

idxs = list(data.obs.index)

# dprint(idxs)

idxs = [x.split('=')[1] for x in idxs]


# dprint(filtered_data)

# dprint(data.var.shape)



# iterate over pucks
pids = [i for i in np.arange(1,208,2)]
pids.remove(5)
pids.remove(77)
pids.remove(167)

pids = pids

out_mat = np.zeros((1,4))

# for idx, pid in enumerate(pids[51:53]):
for idx, pid in enumerate(pids[70:]):
    inp_file = f'{in_folder}/dd{pid}_CTMapping.h5ad'
    dprint(inp_file)
    data = ann.read_h5ad(inp_file)

    # get cell type index
    filtered_idxs = [i for i,x in enumerate(idxs) if 'Nrk' in x]

    # sum up selected rows
    # dprint(data.X.getrow(10))
    selected_data = data.X[filtered_idxs,:]

    aggregated_data = np.sum(selected_data, axis=0)

    higher = (aggregated_data > 0)
    # dprint(higher)
    # dprint(np.sum(higher))

    higher_inds = np.where(higher == True)[1]
    # dprint(higher_inds.shape)

    fvals = aggregated_data[0,higher_inds]
    fvals = np.swapaxes(fvals, 0,1)
    # dprint(fvals.shape)

    # get positions
    # pos_x = data.var
    positions = data.var[['CCF3D_x', 'CCF3D_y', 'CCF3D_z']].values[higher_inds,:]
    # dprint(positions.shape)

    # merge positions and fvals
    merged_mat = np.hstack([positions, fvals])
    # dprint(merged_mat.shape)
    # dprint(merged_mat)

    # add to global writeout matrix
    out_mat = np.vstack((out_mat, merged_mat))


# write output to csv

out_mat = out_mat[1:, :] # remove the initial zero vector
dprint(out_file)
header = 'xx,yy,zz,fval'
np.savetxt(out_file, out_mat, delimiter=',', fmt='%1.3f', header=header, comments='')



# # write to ouptut vtp
# out_file_vtp = f'{out_folder}/fvals.vtp'

# colors = vtk.vtkFloatArray()

# # Create a vtkPoints object and store the points in it
# points = vtk.vtkPoints()
# for i in range(out_mat.shape[0]):
#     points.InsertNextPoint(out_mat[i,0], out_mat[i,1], out_mat[i,2])
#     colors.InsertNextValue(out_mat[i,3])

# # Create a polydata and add the points to it
# polydata = vtk.vtkPolyData()
# polydata.SetPoints(points)
# polydata.GetPointData().SetScalars(colors)


# # Write out the points to a vtp file
# pointsWriter = vtk.vtkXMLPolyDataWriter()
# pointsWriter.SetFileName(out_file_vtp)
# pointsWriter.SetInputData(polydata)
# pointsWriter.Write()

