import os

parent = "/Users/mraj/Desktop/work/data/mouse_atlas/initial_converted_subset"

for i in range(130):
    path =parent+"/"+str(i)
    print path
    try:
        os.mkdir(path)
    except:
        print path+" EXISTS"
