import os
import pandas as pd
import numpy as np
from statistics import mean
import math

def write_backbone_atoms(path_edit, header, atom_bb, ROI_start, ROI_end):
    
    # path_edit = '/Users/adk9hq/Documents/Research_UVA/measure_domain/rasa1'
    pdbfiles_edit = os.listdir(path_edit)
    
    for edit_file in pdbfiles_edit:
        if edit_file.startswith("edit_"):
          
            file_path_edit = f"{path_edit}/{edit_file}"
            N_file = header+edit_file
            newname = os.path.join(path_edit, N_file)   
            N_file_new = open(newname, 'w')

            for line in open(file_path_edit):

                list = line.split()
                atom_id = list[0]
                atom_type = list[1]
                atom = list[2]
                res = list[3]
                chain = list[4]
                resnum = int(list[5])
                x = float(list[6])
                y = float(list[7])
                z = float(list[8])

                if chain == 'A' and atom == atom_bb:
                    if resnum in range(int(ROI_start), int(ROI_end)):
                        N_file_new.write(line)       
           


def cal_distance(path_N_files, header):
    #path_N_files = '/Users/adk9hq/Documents/Research_UVA/measure_domain/rasa1'
    N_files = os.listdir(path_N_files)
    list_files = []

    for i in N_files:
        if i.startswith(header+"edit_"):
            list_files.append(i)

    distance_file = header + '_distance.txt'
    newname = os.path.join(path_N_files, distance_file)   
    f3 = open(newname, 'w')
    f3.write("Distance\tpdb_id\n")

    list_dist = []

    for usefile in list_files:
        f = f"{path_N_files}/{usefile}"

        f1 = open(f, 'r')

        for j in f1:
            list_1 = j.split()
            resnum_1 = int(list_1[5])
            x1 = float(list_1[6])
            y1 = float(list_1[7])
            z1 = float(list_1[8])

            f2 = open(f, 'r')
            for k in f2:
                list_2 = k.split()
                resnum_2 = int(list_2[5])
                x2 = float(list_2[6])
                y2 = float(list_2[7])
                z2 = float(list_2[8])

                distance = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
                list_dist.append(distance)

        PDB_ID = usefile.split('.')[0][-4:]
        
        if list_dist == []:
            print("empty list",PDB_ID)
        else:
            largest_distance = max(list_dist)
            print(header, largest_distance, PDB_ID)   

            f3.write(str(largest_distance)+"\t"+str(PDB_ID)+"\n")
            list_dist.clear()

    f3.close()    

def average_distance(PATH, header):
    
    #PATH = '/Users/adk9hq/Documents/Research_UVA/measure_domain/rasa1'
    distance_file = header + '_distance.txt'
    newname = os.path.join(PATH, distance_file)   

    lines = pd.read_csv(newname, sep="\t")
    mean = lines['Distance'].mean()
    print(header, round(mean,2))