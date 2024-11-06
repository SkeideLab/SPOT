import os
import pandas as pd
import numpy as np
import nibabel as nib
from nilearn import surface
from neuroCombat import neuroCombat

left = pd.read_csv("/data/p_02915/SPOT/covars_hemi-L_tsnr_v1.csv")
left = np.mean(left, axis=1)
right = pd.read_csv("/data/p_02915/SPOT/covars_hemi-R_tsnr_v1.csv")
right = np.mean(right, axis=1)

combined_subject_data = pd.concat([left, right], axis=1)

print(combined_subject_data.shape)
# Create a covariates DataFrame
covars = pd.read_csv("/data/p_02915/SPOT/covars_hemi-L.csv")
save_raw_data = pd.DataFrame(combined_subject_data.T)
save_raw_data.to_csv(f"/data/p_02915/SPOT/raw_tnsr_v1.csv", index=False, header=False)
# Check unique values in covariates
print(covars[['site', 'age', 'sex']].nunique())

# Ensure all necessary columns are of appropriate type
covars['site'] = covars['site'].astype('category')
covars['sex'] = covars['sex'].astype('category')

# Step 4: Apply ComBat harmonization to the correlation data
#combat_harmonized = neuroCombat(dat=combined_subject_data.T,          # Transpose the data
#                                covars=covars[['site', 'age', 'sex']],               # Covariates DataFrame
#                                batch_col='site',            # Adjust for site
#                                continuous_cols=['age'],     # Preserve 'age'
#                                categorical_cols=['sex'])["data"] # Preserve 'gender'
#save_combat_data = pd.DataFrame(combat_harmonized)
#save_combat_data.to_csv(f"/data/p_02915/SPOT/combat_tnsr.csv", index=False, header=False)


