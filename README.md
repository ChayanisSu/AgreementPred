# Description
AgreementPred_Codes folder contains scripts necessary to implement AgreementPred.

# Data
Due to space limitation on Github. The the data necessary to implement AgreementPred should be downloaded through OneDrive 

https://chula-my.sharepoint.com/:f:/g/personal/chayanis_su_chula_ac_th/EuYYZU-8BixIu-Ox3wkLp-gB-s-ezneIsGDkybCHBLQUUw?e=SaO0Br

Transfer all folders into Agreement_Codes main folder prior to implementing the scripts.

The followint data are provided
1. Pubchem annotation of 9721 compounds in Annotated-Compound dataset (Annotated-Compound_Annotations.csv)
2. Side effects of 1376 drugs in Annotated-SE dataset (Annotated-SE_Annotations.csv)
3. ATC classification of 1157 drugs from iSEA (https://pubs.acs.org/doi/10.1021/ci400155x) dataset. ATC class (second-level ATC) is provided in the original publication. The higher-level ATCs were extracted from DrugBank database.
4. Precalculated molecular representations of 3 datasets of 1-3
5. All pubchem compounds used in this study (All_Compound.csv)
6. A molecular representation name list, mapping representation names used in scripts to those used in the article (29FP_lst.csv)
7. Lists of Pubchem compounds contained in each sample dataset used for the study

# Implementation
To perform pharmacological category recommendation using AgreementPred, the following 3 steps must be conducted in a consecutive manner:

1. Molecular representation generation for query compounds
2. MSC search
3. Prediction

# Molecular representation generation
Scripts for molecular representation generation used for AgreementPred are derived from the following sources:
- Boldini D, Ballabio D, Consonni V, Todeschini R, Grisoni F, Sieber SA. Effectiveness of molecular fingerprints for exploring the chemical space of natural products. J Cheminform. 2024;16(1):35. (DOI: https://doi.org/10.1186/s13321-024-00830-3)
- CDK_pywrapper package (https://pypi.org/project/CDK-pywrapper/)
- Torchdrug (https://torchdrug.ai/docs/installation.html)

As each source of molecular representation has some special requirements, the following Jupyter notebook (.ipynb) should be used for molecular representation generation.

1. CDK+RDKIT+Minhash_FP.ipynb
2. JMAP_FP.ipynb
3. InfoGraph.ipynb (contained in Torchdrug_InfoGraph folder)

Each notebook requires different python environment. Please use .yml files provided to create working environment as instructed in the notebooks

# MSC search
Following molecular representation generation, MSC_search.ipynb notebook can be used to perform MSC search for query compounds (Look for "For query compounds" section in the notebook)


# Prediction
Finally, Prediction_Application.ipynb notebook can be used to generate pharmacological category recommendation for each query compound. N of MSCs and threshold parameters can be adjusted as instructed in the scripts.

The prediction results consist of two files:
1. Query_similar_compounds.csv --- specifying MSCs of the query compound calculated using each representation
2. Query_Threshold_XX_Predictions.csv --- containing recommended categories for each query compound. XX in the file name indicates the threshold parameter used for prediction.


# Prediction performance
Scripts and instruction to evaluate prediction performance as performed in the study can be found in Prediction_Performances.ipynb notebook

