### process the original beta files from Loyfer et.al., 2023
### to obtain proportion of methylated reads per CpG site across the genome
### filtered for CpGs with minimum read depth of 10

import numpy as np
import pandas as pd
import os

os.chdir('/path/to/project/folder/')

def processBeta(filename, minRD = 10):
    content = np.fromfile(filename, dtype=np.uint8).reshape((-1, 2))
    dnam = pd.DataFrame(np.where(content[:, 1] > minRD,  content[:, 0]/content[:, 1], np.nan))
    return(dnam)

colanno = pd.read_csv('metadata/colanno27.csv', sep = ',', header=0, index_col=0)       
pathList = colanno.iloc[:,0].values
betas = pd.concat([processBeta(i,10) for i in pathList], axis=1)

betas.to_csv('output/0_loyfer_methylProp27.csv')