#!/usr/bin/env python3

import pandas as pd

# function to read the csv
def read_mlst_csv(csv):
    colnames=['Sample', 'organism', 'ST', 'arcC', 'aroE', 'glpF', 'gmk', 'pta', 'tpy', 'yqiL'] 
    df = pd.read_csv(csv, delimiter= ',', index_col=False, header=None, names=colnames)
    # strip the other text off and remove > and ~
    df['arcC'] = df['arcC'].str.replace("arcC", "").str.strip("()").str.replace("?", "").str.replace("~", "")
    df['aroE'] = df['aroE'].str.replace("aroE", "").str.strip("()").str.replace("?", "").str.replace("~", "")
    df['glpF'] = df['glpF'].str.replace("glpF", "").str.strip("()").str.replace("?", "").str.replace("~", "")
    df['gmk'] = df['gmk'].str.replace("gmk", "").str.strip("()").str.replace("?", "").str.replace("~", "")
    df['pta'] = df['pta'].str.replace("pta", "").str.strip("()").str.replace("?", "").str.replace("~", "")
    df['tpy'] = df['tpy'].str.replace("tpi", "").str.strip("()").str.replace("?", "").str.replace("~", "") # why it is tpi, who knows
    df['yqiL'] = df['yqiL'].str.replace("yqiL", "").str.strip("()").str.replace("?", "").str.replace("~", "")
    # drop organism and ST
    df = df.drop(['organism', 'ST'], axis=1)
    return df




def summarise_contigs(summary_list, output, saureus):
    # read into list       
    summaries = []
    l =summary_list

    for a in l:
        tmp_summary = read_mlst_csv(a)
        summaries.append(tmp_summary)

    # make into combined dataframe
    total_summary_df = pd.concat(summaries,  ignore_index=True)

    # merge s aureus into the summary df
    colnames=['ST', 'arcC', 'aroE', 'glpF', 'gmk', 'pta', 'tpy', 'yqiL', 'clonal_complex'] 
    saureus_df = pd.read_csv(saureus, delimiter= '\t', index_col=False, header=None, names=colnames)
    print(total_summary_df)
    # convert the locus to string for merge 
    saureus_df[['arcC', 'aroE', 'glpF', 'gmk', 'pta', 'tpy', 'yqiL']] = saureus_df[['arcC', 'aroE', 'glpF', 'gmk', 'pta', 'tpy', 'yqiL']].astype(str)
    # convert the total_summary_df 
    total_summary_df['ST']=pd.merge(total_summary_df, saureus_df, how="left", on=['arcC', 'aroE', 'glpF', 'gmk', 'pta', 'tpy', 'yqiL']).ST
    total_summary_df['clonal_complex']=pd.merge(total_summary_df, saureus_df, how="left", on=['arcC', 'aroE', 'glpF', 'gmk', 'pta', 'tpy', 'yqiL']).clonal_complex
    total_summary_df = total_summary_df.fillna('Novel')
    total_summary_df.to_csv(output, sep=",", index=False)

     
summarise_contigs(snakemake.input.mlsts, snakemake.output.out, snakemake.params.saureus)




