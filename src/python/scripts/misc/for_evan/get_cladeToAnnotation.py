"""
Scrpt to generate a map from clade ID to clade annotation


Created by Mukund on 2024-09-27
"""

# directory: mouse_atlas/single_cell/s3/

z = zarr.open('aggedSCdata_240531.zarr')
clades = list(z.metadata.clades)
cladesAnno = list(z.metadata.cladesAnnotations)
zipped = list(zip(clades, cladesAnno))
np.savetxt('cladeToAnnotationMap.csv', zipped, fmt='%s', delimiter=',')
