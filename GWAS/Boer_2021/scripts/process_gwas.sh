

module load r/4.2.2
module load ucsctools
module load python/3.9.6
module load plink


# Create parent directories
for OA in AllOA FingerOA HandOA HipOA KneeHipOA KneeOA SpineOA THR ThumbOA TJR TKR
do
    mkdir -p ${OA}/leads
    mkdir -p ${OA}/summary_stats
done

################################################################ Summary Statistics
for OA in AllOA FingerOA HandOA HipOA KneeHipOA KneeOA SpineOA THR ThumbOA TJR TKR
do

    # Download
    sbatch --wrap="wget -O ${OA}/summary_stats/KP.Format.GO.FILTER.GW.${OA}.FULL.09052019.txt.gz https://personal.broadinstitute.org/ryank/KP.Format.GO.FILTER.GW.${OA}.FULL.09052019.txt.gz"

    # Reformat and split by chromosome
    sbatch -J ${OA}_gwas -t 24:00:00 --mem=8G --wrap="Rscript reformat_summstats.R ${OA}"

    # Lift over

    # Create temporary BED-formatted file from reformatted leads
    for chr in {1..22}
    do
        awk -F ',' '{print $2 "\t" $3 "\t" $3 "\t" $1}' ${OA}/summary_stats/${OA}_chr${chr}.csv > ${OA}/summary_stats/${OA}_chr${chr}_temp.bed

        # Remove header
        sed -i '1d' ${OA}/summary_stats/${OA}_chr${chr}_temp.bed

        # Add 'chr' prefix to chromosome
        awk '$1="chr"$1' ${OA}/summary_stats/${OA}_chr${chr}_temp.bed > ${OA}/summary_stats/${OA}_chr${chr}_temp_chr.bed

        # Add 1 bp for end position
        awk '{$3+=1}1' ${OA}/summary_stats/${OA}_chr${chr}_temp_chr.bed > ${OA}/summary_stats/${OA}_chr${chr}_temp_end.bed

        # Lift over to hg38
        sbatch -J ${OA}_liftover --wrap="liftOver -bedPlus=1 ${OA}/summary_stats/${OA}_chr${chr}_temp_end.bed \
            /proj/phanstiel_lab/References/genomes/hg19/hg19ToHg38.over.chain.gz ${OA}/summary_stats/${OA}_chr${chr}_liftOver.bed \
            ${OA}/summary_stats/${OA}_chr${chr}_liftOver_unmapped.bed" 

        # Join back with original summary stats
        sbatch -J ${chr}_${OA}_join -t 12:00:00 --mem=8G --wrap="Rscript join_summstats_liftOver.R ${OA} ${chr}"

        # Get rsIDs(?)
    done

done

######################################### Leads

# 1) Parse leads from supp table, separating the leads by OA type and getting the relevant info/ids/hg19position
sbatch -t 4:00:00 -J parseLeads --wrap="Rscript parseLeads.R Boer_gwas_sup3.csv"

for OA in AllOA FingerOA HandOA HipOA KneeHipOA KneeOA SpineOA THR ThumbOA TJR TKR
do
    # 2) Get most updated rsID for lead and get proper alleles for D/I
    sbatch -t 12:00:00 --mem=16G -J rsID_leads --wrap="python3 lead_rsids.py ${OA}/leads/${OA}_leads.csv ${OA} /work/users/n/e/nekramer/References/GRCh37.p13/dbSNP dbSNP155.GRCh37.p13"

    # 3) Lift Over
    
    # Create temporary BED-formatted file from reformatted leads
    awk -F ',' '{print $2 "\t" $3 "\t" $3 "\t" $1}' ${OA}/leads/${OA}_leads_rsid.csv > ${OA}/leads/${OA}_leads_temp.bed

    # Remove header
    sed -i '1d' ${OA}/leads/${OA}_leads_temp.bed

    # Add 'chr' prefix to chromosome
    awk '$1="chr"$1' ${OA}/leads/${OA}_leads_temp.bed > ${OA}/leads/${OA}_leads_temp_chr.bed

    # Add 1 bp for end position
    awk '{$3+=1}1' ${OA}/leads/${OA}_leads_temp_chr.bed > ${OA}/leads/${OA}_leads_temp_end.bed

    # Lift over to hg38
    sbatch -J ${OA}_liftover --wrap="liftOver -bedPlus=1 ${OA}/leads/${OA}_leads_temp_end.bed \
        /proj/phanstiel_lab/References/genomes/hg19/hg19ToHg38.over.chain.gz ${OA}/leads/${OA}_leads_liftOver.bed \
         ${OA}/leads/${OA}_leads_liftOver_unmapped.bed"

    # 4) Join lifted over leads back with original data
    sbatch -t 6:00:00 -J ${OA}_joinLiftOver --wrap="Rscript join_liftOver.R ${OA} ${OA}/leads/${OA}_leads_liftOver.bed"

    # 5) Try getting hg38 position based on rsID for any variants that couldn't be lifted over
    sbatch -t 12:00:00 --mem=16G -J leads_hg38 --wrap="python3 rsid_hg38.py ${OA}/leads/${OA}_leads_liftOver.csv ${OA} proj/phanstiel_lab/References/genomes/GENCODE.GRCh38.p13/dbSNP dbSNP155.GRCh38.p13"

    # 6) Make directories for LD buddies
    mkdir -p ${OA}/leads/ALL_ld
    mkdir -p ${OA}/leads/EUR_ld
    mkdir -p ${OA}/leads/ALL_1000G_snps
    mkdir -p ${OA}/leads/EUR_1000G_snps

    # 7) For each lead SNP, get LD buddies from 1000G (using chrom:hg19POS to match 1000G IDs)
    for snp in `awk -F ',' '{print $5}' ${OA}/leads/${OA}_leads_liftOver_final.csv | awk 'NR!=1 {{print}}'`
    do
        # Get chr
        chr=$(echo $snp | awk -F ':' '{print $1}')    
        # Get ALL LD buddies
        sbatch -t 12:00:00 --mem=16G --wrap="plink --bfile /proj/phanstiel_lab/References/genomes/1000G/GRCh37/ALL/all_snps/1000G_phase3_chr${chr}_noid --ld-snp ${snp} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out ${OA}/leads/ALL_ld/${snp}_ld"
        # Get EUR LD buddies
        sbatch -t 12:00:00 --mem=16G --wrap="plink --bfile /proj/phanstiel_lab/References/genomes/1000G/GRCh37/EUR/all_snps/EUR_1000G.GRCh37_chr${chr}_noid --ld-snp ${snp} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out ${OA}/leads/EUR_ld/${snp}_ld"

        # Get alleles for LD buddies from 1000G
        awk -F ' ' '{print $6}' ${OA}/leads/ALL_ld/${snp}_ld.ld | awk 'NR!=1 {{print}}' > ${OA}/leads/ALL_ld/${snp}_ld_buddies.txt
        sbatch -t 12:00:00 --mem=8G --wrap="plink --bfile /proj/phanstiel_lab/References/genomes/1000G/GRCh37/ALL/all_snps/1000G_phase3_chr${chr}_noid --extract ${OA}/leads/ALL_ld/${snp}_ld_buddies.txt --make-bed --out ${OA}/leads/ALL_1000G_snps/${snp}_1000G"
        
        awk -F ' ' '{print $6}' ${OA}/leads/EUR_ld/${snp}_ld.ld | awk 'NR!=1 {{print}}' > ${OA}/leads/EUR_ld/${snp}_ld_buddies.txt
        sbatch -t 12:00:00 --mem=8G --wrap="plink --bfile /proj/phanstiel_lab/References/genomes/1000G/GRCh37/EUR/all_snps/EUR_1000G.GRCh37_chr${chr}_noid --extract ${OA}/leads/EUR_ld/${snp}_ld_buddies.txt --make-bed --out ${OA}/leads/EUR_1000G_snps/${snp}_1000G"

        # Join LD buddies with lead data and get alleles for buddies
        sbatch -t 12:00:00 --mem=4G --wrap="Rscript reformat_LDbuddies.R ${snp} ${OA} ALL ${OA}/leads/ALL_ld/${snp}_ld.ld"
        sbatch -t 12:00:00 --mem=4G --wrap="Rscript reformat_LDbuddies.R ${snp} ${OA} EUR ${OA}/leads/EUR_ld/${snp}_ld.ld"

    done

    # 8) Concatenate all leads_ldbuddies, making sure to keep leads that weren't found in 1000G
    sbatch --wrap="Rscript concatenate_lead_lds.R ${OA} ALL"
    sbatch --wrap="Rscript concatenate_lead_lds.R ${OA} EUR"

    # 9) Get rsIDs for LD buddies from hg19 position
    sbatch -t 72:00:00 --mem=32G -J ${OA}_ALL_rsid_ldbuddies --wrap="python3 ld_rsids.py ${OA}/leads/ALL_HipOA_leads_ld.csv ${OA} ALL /work/users/n/e/nekramer/References/GRCh37.p13/dbSNP dbSNP155.GRCh37.p13"
    sbatch -t 72:00:00 --mem=32G -J ${OA}_EUR_rsid_ldbuddies --wrap="python3 ld_rsids.py ${OA}/leads/EUR_HipOA_leads_ld.csv ${OA} EUR /work/users/n/e/nekramer/References/GRCh37.p13/dbSNP dbSNP155.GRCh37.p13"
    
    # 10) Lift over LD buddies

    # Create temporary bed-formatted file parsing LD buddy IDs
    awk -F ',' '{split($51,a, ":"); print a[1] "\t" a[2] "\t" a[2] "\t" $53}' ${OA}/leads/ALL_${OA}_leads_ld_rsid.csv > ${OA}/leads/ALL_${OA}_leads_ld_temp.bed
    awk -F ',' '{split($51,a, ":"); print a[1] "\t" a[2] "\t" a[2] "\t" $53}' ${OA}/leads/EUR_${OA}_leads_ld_rsid.csv > ${OA}/leads/EUR_${OA}_leads_ld_temp.bed

    # Remove header
    sed -i '1d' ${OA}/leads/ALL_${OA}_leads_ld_temp.bed
    sed -i '1d' ${OA}/leads/EUR_${OA}_leads_ld_temp.bed

    # Add 'chr' prefix to chromosome
    awk '$1="chr"$1' ${OA}/leads/ALL_${OA}_leads_ld_temp.bed > ${OA}/leads/ALL_${OA}_leads_ld_temp_chr.bed
    awk '$1="chr"$1' ${OA}/leads/EUR_${OA}_leads_ld_temp.bed > ${OA}/leads/EUR_${OA}_leads_ld_temp_chr.bed

    # Add 1 bp for end position
    awk '{$3+=1}1' ${OA}/leads/ALL_${OA}_leads_ld_temp_chr.bed > ${OA}/leads/ALL_${OA}_leads_ld_temp_end.bed
    awk '{$3+=1}1' ${OA}/leads/EUR_${OA}_leads_ld_temp_chr.bed > ${OA}/leads/EUR_${OA}_leads_ld_temp_end.bed

    # Lift over to hg38
    sbatch --wrap="liftOver -bedPlus=1 ${OA}/leads/ALL_${OA}_leads_ld_temp_end.bed \
        /proj/phanstiel_lab/References/genomes/hg19/hg19ToHg38.over.chain.gz ${OA}/leads/ALL_${OA}_leads_ld_liftOver.bed \
         ${OA}/leads/ALL_${OA}_leads_ld_liftOver_unmapped.bed"

    sbatch --wrap="liftOver -bedPlus=1 ${OA}/leads/EUR_${OA}_leads_ld_temp_end.bed \
        /proj/phanstiel_lab/References/genomes/hg19/hg19ToHg38.over.chain.gz ${OA}/leads/EUR_${OA}_leads_ld_liftOver.bed \
         ${OA}/leads/EUR_${OA}_leads_ld_liftOver_unmapped.bed"

    # Join back with original lead ld data
    sbatch --wrap="Rscript join_liftOver_leadld.R ${OA} ALL ${OA}/leads/ALL_${OA}_leads_ld_liftOver.bed"
    sbatch --wrap="Rscript join_liftOver_leadld.R ${OA} EUR ${OA}/leads/EUR_${OA}_leads_ld_liftOver.bed"

    # 11) Get summary statistic p-values for LD buddies for OA type
    sbatch -t 72:00:00 --mem=8G -J ${OA}_ALL_sum --wrap="Rscript ldbuddy_summarystats.R ${OA} ALL"
    sbatch -t 72:00:00 --mem=8G -J ${OA}_EUR_sum --wrap="Rscript ldbuddy_summarystats.R ${OA} EUR"

done


