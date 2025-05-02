import os
import sys

import pandas as pd
import numpy as np
import json

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools

from CDK_pywrapper import CDK, FPType
from FP_calc.rdkit_fps import *
from FP_calc.minhash_fps import *
from FP_calc.cdk_fps import *

from tqdm import tqdm

def Calc_fps(compound_df, main_folder):

    all_chem = compound_df

    # Creating molecular object from Smiles
    PandasTools.AddMoleculeColumnToFrame(all_chem,'canonicalsmiles','mol')
    df = all_chem.reset_index(drop=True)
    mols = df['mol'].tolist()
    smiles = list(df['canonicalsmiles'])

    #CDKFP
    cdk = CDK(fingerprint=FPType.FP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/cdkFP.csv")

    #ExtFP
    cdk = CDK(fingerprint=FPType.ExtFP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/ExtFP.csv")

    #EStateFP
    cdk = CDK(fingerprint=FPType.EStateFP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/EStateFP.csv")

    #GraphFP
    cdk = CDK(fingerprint=FPType.GraphFP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/GraphFP.csv")

    #MACCSFP
    cdk = CDK(fingerprint=FPType.MACCSFP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/maccsFP.csv")

    #PubchemFP
    cdk = CDK(fingerprint=FPType.PubchemFP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/PubchemFP.csv")

    #KRFP
    cdk = CDK(fingerprint=FPType.KRFP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/KRFP.csv")

    #AP2DFP
    cdk = CDK(fingerprint=FPType.AP2DFP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/AP2DFP.csv")

    #HybridFP
    cdk = CDK(fingerprint=FPType.HybridFP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/HybridFP.csv")

    #LingoFP
    cdk = CDK(fingerprint=FPType.LingoFP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/LingoFP.csv")

    #SPFP
    cdk = CDK(fingerprint=FPType.SPFP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/SPFP.csv")

    #CircFP
    cdk = CDK(fingerprint=FPType.CircFP)
    molfp = cdk.calculate(mols)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/CircFP.csv")

    #Daylight
    fp = calc_DAYLIGHT(smiles)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/cdkDaylight.csv")

    #RDKIT
    #EC1024
    fp = calc_ECFP(mols)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/rdkitECFP.csv")

    #EC2048
    fp = calc_Morgan(mols)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/MorganFP.csv")

    #FC1024
    fp = calc_FCFP(mols)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/rdkitFCFP.csv")
    
    #Avalon
    fp = calc_AVALON(mols)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/rdkitAvalon.csv")

    #TT
    fp = calc_TT(mols)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/rdkitTT.csv")

    #AP
    fp = calc_AP(mols)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/rdkitAP.csv")

    #RDKIT
    fp = calc_RDKIT(mols)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/rdkitRDKIT.csv")

    #Minhash
    #MAP4
    fp = calc_MAP4(mols)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/minhashMAP4.csv")

    #MHFP
    fp = calc_MHFP(mols)
    molfp = pd.DataFrame(fp)
    molfp['cid'] = df['cid']
    molfp['cmpdname'] = df['cmpdname']
    molfp.to_csv(f"{main_folder}/Representations/minhashMHFP.csv")

