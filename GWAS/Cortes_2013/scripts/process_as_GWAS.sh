module load r/4.3.1
module load ucsctools
module load samtools
module load plink


####################### Harmonized summary stats from the GWAS catalog

# Harmonized files are already in GRCh38 genome build 

sbatch --wrap="wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005529/harmonised/23749187-GCST005529-EFO_0003898.h.tsv.gz"


################### Leads taken from Table 1 (non-MHC loci) of Cortes et al. 2013

# positions are in NCBI build 36 -> lift over 

# 1) lift over to hg38
# 1a) Create temporary BED-formatted file from leads file
awk -F ',' '{print $2 "\t" $3 "\t" $3 "\t" $1}' IGAS_2013_Table1.csv > IGAS_leads_temp.bed
# 1b) Remove header
sed -i '1d' IGAS_leads_temp.bed
# 1c) Add 1 bp for end position
awk '{$3+=1}1' IGAS_leads_temp.bed > IGAS_leads_temp_end.bed

# 5aiv) Lift over to hg38
sbatch --wrap="liftOver -bedPlus=1 IGAS_leads_temp_end.bed \
    /proj/phanstiel_lab/References/genomes/hg38/hg18ToHg38.over.chain.gz IGAS_leads_liftOver.bed \
    IGAS_leads_liftOver_unmapped.bed"

# Rejoin lift over to original lead data
sbatch --wrap="Rscript join_liftOver.R"

#2) get LD buddies for leads 

# Trying one verison of variant ID
for snp in `awk -F ',' '{print $2}' IGAS_2013_leads_hg38.csv | awk 'NR!=1 {{print}}'`
do
    # Get chrom
    chr=$(echo $snp | awk -F ':' '{print $1}')

    # get ld from 1000G EUR hg38
    sbatch --mem=4G -t 12:00:00 --wrap="plink --bfile /proj/phanstiel_lab/References/genomes/1000G/GRCh38/EUR/all_snps/EUR_1000G.GRCh38.20181129_chr${chr} --ld-snp ${snp} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out ${snp}_ld"

done

# Trying second version of variant ID
for snp in `awk -F ',' '{print $3}' IGAS_2013_leads_hg38.csv | awk 'NR!=1 {{print}}'`
do
    # Get chrom
    chr=$(echo $snp | awk -F ':' '{print $1}')

    # get ld from 1000G EUR hg38
    sbatch --mem=4G -t 12:00:00 --wrap="plink --bfile /proj/phanstiel_lab/References/genomes/1000G/GRCh38/EUR/all_snps/EUR_1000G.GRCh38.20181129_chr${chr} --ld-snp ${snp} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out ${snp}_ld"

done

# Join in and reformat LD buddies and summary stats
sbatch -t 12:00:00 --mem=16G --wrap="Rscript join_reormat_ld_buddies.R"