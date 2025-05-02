import os
import sys

import pandas as pd
import numpy as np
import json

from FP_calc.jmap_fps import *

from tqdm import tqdm

def Calc_JMAP(compound_df, main_folder):

    all_chem = compound_df

    # Creating molecular object from Smiles
    df = all_chem.reset_index(drop=True)
    smiles = list(df['canonicalsmiles'])

    #DFS
    fp = calc_DFS(smiles)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/jmapDFS.csv")

    #ASP
    fp = calc_ASP(smiles)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/jmapASP.csv")

    #LSTAR
    fp = calc_LSTAR(smiles)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/jmapLSTAR.csv")

    #RAD2D
    fp = calc_RAD2D(smiles)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/jmapRAD2D.csv")

    #PH2
    fp = calc_PH2(smiles)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/jmapPH2.csv")

    #PH3
    fp = calc_PH3(smiles)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/jmapPH3.csv")

