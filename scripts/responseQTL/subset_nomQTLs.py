import pandas as pd

with open('output/reQTL/chr16nominal_results.csv', "a+") as output:
  # Write header
    output.write('gene_id gene_chr gene_start gene_end gene_strand num_cis_variants dist variantID variant_chr variant_start variant_end nom_pval r_squared beta best_hit pval_nominal_threshold\n')
    with open('/pine/scr/n/e/nekramer/iteratePEER_RNAKitBatch_DNAKitBatch/output/qtl/FNF_PEER_k10_genoPC_RNAKitBatch_DNAKitBatch_nom1Mb.txt', "r") as f:
      for line in f:
        qtldata = line.rstrip().split(" ")
        variant_chr = qtldata[8]
        variant_pos = int(qtldata[9])
        
        if variant_chr == 'chr16' and (variant_pos >= 69859821 or variant_pos <= 70736995):
          output.write(line.rstrip() + "\n")
