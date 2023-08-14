import pandas as pd
from utils import getRsids
import sys


reqtls = sys.argv[1]

reqtl_positions = pd.read_csv(reqtls)
reqtl_rsIDs = getRsids(reqtl_positions)


reqtl_positions_rsIDs = pd.concat([reqtl_positions, reqtl_rsIDs], axis = 1)
reqtl_positions_rsIDs.to_csv('output/reQTL/FNF_sig_reQTLs_rsids.csv', index = False)
