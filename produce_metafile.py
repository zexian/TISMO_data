import os
import pandas as pd

def main():
    files_df = pd.DataFrame(columns = ["SampleName", "GSE_ID", "Cell_Status", "Sequencing_Type", "File_Path"])
    annot_df = pd.read_csv("may_vivo2_srx.csv") #input csv with 2 columns: GSE_ID, SRX_ID
    path = '/liulab/cheryl/data/may21_icb/vivo/'
    gse_list = []
    for direct in os.listdir(path):
        if 'GSE' in direct:
            gse_list.append(direct)
    #print(len(gse_list))        
    
    i = 0
    num_dups = 0
    dup_list = []
    not_found = []
    for gse in gse_list:
        print(sorted(os.listdir(os.path.join(path, gse))))
        for srx in sorted(os.listdir(os.path.join(path, gse))): 
            print("SRX:", srx)
            
            srx_loc = annot_df.SRX_ID[annot_df.SRX_ID == srx].index.to_list()
            print("SRX_LOC:",srx_loc)
            
            '''if len(srx_loc) > 1:
                #num_dups +=1
                dup_list.append(srx)'''
            index = None
            if srx_loc:
                files_df.loc[i, "SRX_ID"] = srx
                for ind in srx_loc:
                    if annot_df.loc[ind,"GSE_ID"] == gse:
                        index = ind
                #index = srx_loc[0]
                gse_id = annot_df.loc[index,"GSE_ID"]
                #print(gse)
                gse_id = gse_id.strip()
                files_df.loc[i, "GSE_ID"] = gse_id
                files_df.loc[i, "SampleName"] = srx
                if 'vivo' in path:
                    files_df.loc[i, "Cell_Status"] = "in_vivo"
                elif 'vitro' in path:
                    files_df.loc[i, "Cell_Status"] = "in_vitro"   
                #print(cell_status, sample_name)
            

            else:
                print(srx,"not found")
                not_found.append(srx)
                
            srx_path = os.path.join(os.path.join(path, gse),srx)
            if srx_loc:
                for srr in os.listdir(srx_path):
                    fastq_list = os.listdir(os.path.join(srx_path,srr))
                    #print(len(fastq_list), fastq_list[0])
                    #three_options = 0
                    beg_path = os.path.join(srx_path,srr) + '/'
                    if len(fastq_list) == 3:
                        #print("here")
                        files_df.loc[i, "Sequencing_Type"] = "Paired"
                        pe_list = []
                        for f in fastq_list:
                            if '_' in f:
                                pe_list.append(f)
                        fastq_path = beg_path + pe_list[0] + ';' + beg_path + pe_list[1] 
                        files_df.loc[i, "File_Path"] = fastq_path       

                    elif len(fastq_list) == 2:
                        files_df.loc[i, "Sequencing_Type"] = "Paired"
                        fastq_path = beg_path + fastq_list[0] + ';' + beg_path + fastq_list[1]
                        files_df.loc[i, "File_Path"] = fastq_path
                    elif len(fastq_list) == 1:
                        files_df.loc[i, "Sequencing_Type"] = "Single"
                        fastq_path = beg_path + fastq_list[0]
                        files_df.loc[i, "File_Path"] = fastq_path

                    else:
                        files_df.loc[i, "File_Path"] = "None"  
                i += 1            
    print(files_df.shape[0])
    # output meta.csv
    files_df.to_csv("may_vivo2_meta.csv", encoding='utf-8', index=False)
    print(len(not_found))

if __name__== "__main__":
    main()   

                        
                                


