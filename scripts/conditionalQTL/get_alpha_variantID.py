import sys

with open(sys.argv[2], "a+") as output:
    # Write header
    output.write('rsID,variantID,gene_id,gene_symbol,gene_chr,gene_start,gene_end,gene_strand,variant_chr,variant_start,variant_end,signal,backward_nom_pvalue,backward_r_squared,backward_slope,backward_best_hit,ld_variantID,R2,ld_rsID,variantID_alpha,ld_variantID_alpha\n')
    with open(sys.argv[1],"r") as file:
        for line in file:
            # Skip header
            if line[:4] == "rsID": continue

            qtl_data = line.rstrip().split(",")
            # Convert variantID to alphabetical allele order
            chr,pos,a1,a2 = qtl_data[1].split(":")
            alpha_variantID = "%s:%s:%s:%s" % (chr,pos,sorted([a1,a2])[0],sorted([a1,a2])[1])

            # Convert ld_variantID to alphabetical allele order
            if len(qtl_data) == 19:
                ld_chr,ld_pos,ld_a1,ld_a2 = qtl_data[16].split(":")
                alpha_ldvariantID = "%s:%s:%s:%s" % (ld_chr,ld_pos,sorted([ld_a1,ld_a2])[0],sorted([ld_a1,ld_a2])[1])
                new_line = ",".join(qtl_data) + "," + alpha_variantID + "," + alpha_ldvariantID + "\n"
            elif len(qtl_data) == 20:
                # This probably has some split gene name so we need to join it back together
                gene_symbol = qtl_data[3] + "_" + qtl_data[4]
                ld_chr,ld_pos,ld_a1,ld_a2 = qtl_data[17].split(":")
                alpha_ldvariantID = "%s:%s:%s:%s" % (ld_chr,ld_pos,sorted([ld_a1,ld_a2])[0],sorted([ld_a1,ld_a2])[1])
                new_line = ",".join(qtl_data[0:3]) + "," + gene_symbol + "," + ",".join(qtl_data[5:20]) + "," + alpha_variantID + "," + alpha_ldvariantID + "\n"

            # Make new file line
            
            output.write(new_line)
