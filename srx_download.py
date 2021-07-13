import pandas as pd
import subprocess 
import os
import sys

# ensure enaBrowserTools is downloaded from https://github.com/enasequence/enaBrowserTools
sys.path.insert(1, "/liulab/cheryl/enaBrowserTools/python3")

def main():
    base_path = "/liulab/cheryl/data/may21_icb/vivo"
    dest_dir = ''
    annot_table = pd.read_csv("may_vivo2_srx.csv") # csv with 2 columns: GSE_ID, SRX_ID
    for i in range(len(annot_table)):
        new_gse_num = annot_table.loc[i, 'GSE_ID'].strip()
        print(new_gse_num)
        if new_gse_num != dest_dir: 
            dest_dir = new_gse_num
            if os.path.isdir(dest_dir) == False:
                os.mkdir(dest_dir)
            #num_studies += 1
        srx_id = annot_table.loc[i, 'SRX_ID'].strip()
        print(srx_id)
        full_path = os.path.join(base_path, dest_dir) 
        print(os.path.join(full_path, srx_id))
        # download fastq using enaBrowserTools
        os.system('python /liulab/cheryl/enaBrowserTools/python3/enaDataGet.py -f fastq -d ' + full_path + ' ' + srx_id)
    return   

if __name__== "__main__":
    main()   
