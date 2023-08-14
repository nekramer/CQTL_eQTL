import pandas as pd
from utils import getRsids


variantPositions = pd.read_csv("/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/variants_opposite_effects.csv")
variantRsids = getRsids(variantPositions)

variantsAll = pd.concat([variantPositions, variantRsids], axis = 1)
variantsAll.to_csv('/work/users/n/e/nekramer/Data/CQTL/eQTL/freeze/qtl/variants_opposite_effects_rsids.csv', index = False)


