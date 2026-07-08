# process the bta files

import numpy as np
import pandas as pd
import os

# os.chdir('/SAN/ghlab/epigen/Hiruni/')
# os.chdir('/Users/hirunidb/Library/CloudStorage/OneDrive-Personal/Projects/Hellenthal-SilverLab_LIDo/hellenthalSilverLab_DNAm_pvt/')
os.chdir('/Users/hirunidb/Library/CloudStorage/OneDrive-Personal/Projects/Hellenthal-SilverLab_LIDo/hellenthalSilverLab_DNAm_pvt/loyfer_pancreas/')

def processBeta(filename, minRD = 10):
    content = np.fromfile(filename, dtype=np.uint8).reshape((-1, 2))
    dnam = pd.DataFrame(np.where(content[:, 1] > minRD,  content[:, 0]/content[:, 1], np.nan))
    return(dnam)

colanno = pd.read_csv('metadata/colanno.csv', sep = ',', header=0, index_col=0)       
pathList = colanno.iloc[:,0].values
betas = pd.concat([processBeta(i,10) for i in pathList], axis=1)

betas.to_csv('output/0_loyfer_methylProp36.csv') # blood 
betas.to_csv('output/0_loyfer_methylProp27.csv') # blood 
betas.to_csv('output/0_loyfer_methylProp.csv') # pancreas
