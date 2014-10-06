import pandas as pd
import os

df1 = '/media/toomas/DATA1/aligment_bowtie/MazF_2h/aligned/MazF2h_PNK_index.csv'
df2 = '/media/toomas/DATA1/aligment_bowtie/MazF_2h/aligned/MazF2h_none_index.csv'

def pandas_processer(index):

    abspath = os.path.abspath(index)
    dir_name = os.path.dirname(abspath)
    os.chdir(dir_name)

    #Get root filename for CSV file.
    base_name = os.path.basename(index)
    #Get root and extention of the reference genome fasta file.
    root_ext = os.path.splitext(base_name)
    #Get the base of input file
    genome_index = root_ext[0]
    #Get the file extention input file
    name = genome_index + "_7operons.csv"
    #Creates dataframe from csv
    df = pd.read_csv(index)
    
    #Fills empty dataspaces with 0-s
    df = df.fillna(0)

    #Create a new column in dataframe, which contains the number of rRNA operons which read maps to. 
    df['gene_count'] = df['rrsA'] + df['rrsB'] + df['rrsC'] + df['rrsD'] + df['rrsE'] + df['rrsG'] + df['rrsH']
    #Creates a new fataframe which contains only reads which map to only a certain numbero f rRNA operons
    #selected_df = df.query('gene_count == 7')
    for k in ['rrsA', 'rrsB', 'rrsC', 'rrsD', 'rrsE', 'rrsG', 'rrsH']:
        del df[k]
        
    selected_df = df[df['gene_count'] == 7]
    #Save the new dataframe
    selected_df.to_csv(name,index=True,header=True)


pandas_processer(df1)