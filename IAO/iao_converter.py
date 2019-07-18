"""Nicolas Gampierakis"""

import pandas as pd
import numpy as np
import re

# Read in raster files

# midraster = pd.DataFrame
path_file = "./data/mid/mid_raster_y2d.txt"
# midraster = pd.read_csv(path_file, header=None, sep="\t")
# print(midraster)
# f = open(path_file, 'r+b')
# # read entire content of file into memory
# f_content = f.read()
# f.seek(0)
#     # basically match middle line and replace it with itself and the extra line
# z = str
# z = re.sub(r"[ ,.\t]", r"\t", z)
# f.seek(0)
# # clear file content
# f.truncate()
# # re-write the content with the updated content
# f.write(f_content)
# # close file
# f.close()


output = open("./data/output.txt","w")
input = open(path_file)

for line in input:
    output.write(re.sub(r'\n ', r'\n', line))
    output.write(re.sub(r'  ', r'\t', line))

input.close()
output.close()