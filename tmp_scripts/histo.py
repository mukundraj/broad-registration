"""

python tmp_scripts/histo.py \
"""
import numpy as np
import matplotlib.pyplot as plt


file1 = '/Users/mraj/Downloads/test_data2_gene_jsons_puck189_gene_Gad2.csv'
file2 = '/Users/mraj/Downloads/test_data2_gene_jsons_puck189_gene_Pcp4.csv'

x = np.log(np.genfromtxt(file1, skip_header=1)+1)
x2 = np.log(np.genfromtxt(file2, skip_header=1)+1)

print(np.shape(x))

plt.yscale('log', nonpositive='clip')
plt.hist(x)
plt.show()
plt.yscale('log', nonpositive='clip')
plt.hist(x2)

plt.show()
