"""
Misc helper functions

Created by Mukund on 2022-11-11

"""

"""
Get the gets the nissl index after correcting for offset in pids [175 - 183]

Created by Mukund of 2022-11-11

"""
def get_ofixed_nis_idx(nis_idx):
    if (nis_idx==175):
        corrected_nis_idx = 173
    elif (nis_idx==177):
        corrected_nis_idx = 175
    elif (nis_idx==179):
        corrected_nis_idx = 177
    elif (nis_idx==181):
        corrected_nis_idx = 179
    elif (nis_idx==183):
        corrected_nis_idx = 181
    else:
        corrected_nis_idx = nis_idx
    return corrected_nis_idx


