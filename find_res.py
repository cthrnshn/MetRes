import urllib.request
import os
import biopandas.pdb as bpd
import pandas as pd
import numpy as np
import csv

# creates a new folder to place pdb files if it does not already exist #
def check_directory():
    if not os.path.exists("{}pdb_files".format(newpath)):
        os.makedirs("{}pdb_files".format(newpath))
    return
# from the input fasta file, extracts pdb ids into list #
def sort_fasta(fName):
    all_ids = []
    with open(fName) as file:
        for line in file:
            if line.startswith(">") and line[1:5] not in all_ids:
                all_ids.append(line[1:5])
    all_ids.sort()
    return all_ids
# downloads all pdb files from all_ids list if available, bypasses (+ records id) if not available #
def download_pdb_files(val):
    global omitted_pdb_ids
    global downloaded
    downloaded = True
    if os.path.exists("{}pdb_files\\{}.pdb".format(newpath,val)) == False:
        try:
            url = "https://files.rcsb.org/download/{}.pdb".format(val)
            path = urllib.request.urlretrieve(url, "{}pdb_files\\{}.pdb".format(newpath,val))
        except (urllib.error.HTTPError):
            omitted_pdb_ids.append(val)
            print("OMITTED: " + val)
            downloaded = False
    '''with open("omitted_pdb_ids.csv","w" ) as output:
        for item in omitted_pdb_ids:
            output.write(str(item) + "\n")'''
    return 
# initializes all atoms with record type HETATM and element symbol "m" into a single dataframe, keeps only relevant columns #
def initialize_metals(m):
    m_df = pd.DataFrame()
    df = ppdb.df['HETATM'][(ppdb.df['HETATM']['element_symbol'] == '{}'.format(m))]
    df = df.filter(['atom_number','atom_name','residue_name','residue_number','chain_id','x_coord', 'y_coord', 'z_coord'])
    df.insert(0, "PDB_ID", row, True)
    df = df.set_index("PDB_ID")
    return df
# initializes all atoms with record type ATOM into a single dataframe, keeps only relevant columns #
def initialize_res():
    r_df = pd.DataFrame()
    df = ppdb.df['ATOM']
    df = df.filter(['atom_number','atom_name','residue_name','residue_number','chain_id','x_coord', 'y_coord', 'z_coord'])
    df.insert(0, "PDB_ID", row, True)
    r_df = r_df.append(df, ignore_index = True)
    r_df = r_df.set_index("PDB_ID")
    return r_df
# filters all atoms that fit a specified condition from a given array #
def res_of_interest(cond, all_res):
    idx = np.flatnonzero(cond)
    lst = idx.tolist()
    res = all_res.iloc[lst]
    return res
    
## TOP LEVEL ##
metal_df = pd.DataFrame()
res_df = pd.DataFrame()
all_close_res = pd.DataFrame()
close_res = pd.DataFrame()
omitted_pdb_ids = []
downloaded = True
all_uniq_cool_res = pd.DataFrame()
uniq_cool_res_list = []

## MAIN LEVEL ##
'''
fasta_file = input("Type name of fasta file (include .fasta): ")
newpath = input("Type path of fasta file: ")
r = float(input("Set r: "))
m = input("Type metal element symbol: ").upper()
print(m)
res = input("Type residue of interest: ").upper()
'''
fasta_file = "nr_cu_syn_98percent.fasta"
newpath = ".\Cu\98Percent\\"
r = float("4.0")
m = "CU"
res = "GLN"
atom = "SD"

check_directory()
ids = sort_fasta("{}{}".format(newpath,fasta_file))
count = 0
for row in ids:
    print(row)
    count = count + 1
    print (count)
    download_pdb_files(row)
    if downloaded == True:
        ppdb = bpd.PandasPdb().read_pdb("{}pdb_files\\{}.pdb".format(newpath, row))
        m_df = initialize_metals(m)
        r_df = initialize_res()
        # calculates r of all atoms from all metals #
        for i in range(len(m_df)):
            metal_coords = np.array(m_df.iloc[i,5:8])
            metal_coords = metal_coords.astype(np.float64)
            res_coords = np.array(r_df[['x_coord', 'y_coord', 'z_coord']])
            sub_coords = res_coords - metal_coords
            dist = np.linalg.norm(sub_coords, axis=1)
            r_df['dist_from_metal'] = dist
            close_res = res_of_interest(dist <= r, r_df)
            all_close_res = all_close_res.append(close_res, ignore_index = False)
all_close_res.to_csv(r'.\\res_output.csv')
# filters for only cool_res #
all_close_cool_res = res_of_interest(all_close_res['residue_name'] == "{}".format(res), all_close_res)
cool_res_fName = ".\\{}_output.csv".format(res.lower())
all_close_cool_res.to_csv(r'.\\{}_output.csv'.format(res.lower()))

'''
print("\nTHESE ARE ALL THE CLOSE RESIDUES:")
print(all_close_res)
print("\nTHESE ARE ALL THE CLOSE {} RESIDUES:".format(res))
print(all_close_cool_res)


# creates a csv for all unique_cool_res #
for i in range(len(all_close_cool_res)):
    uniq_cool_res = all_close_cool_res.iloc[[i],[1,2,3,4]]
    data = str(uniq_cool_res)
    if atom in data and data not in uniq_cool_res_list:
        uniq_cool_res_list.append(str(uniq_cool_res))
        uniq_cool_res = uniq_cool_res.filter(['residue_name','residue_number','chain_id'])
        all_uniq_cool_res = all_uniq_cool_res.append(uniq_cool_res, ignore_index = False)
all_uniq_cool_res.to_csv(r'.\\uniq_{}_list.csv'.format(res.lower()))
'''

with open("omitted_pdb_ids.csv","w" ) as output:
    for item in omitted_pdb_ids:
        output.write(str(item) + "\n")
print("Done!")     
