module load r/4.3.1
module load ucsctools
module load samtools
module load plink

# Create parent directories

for RA in AllRA SeropositiveRA
do
    for ancestry in multi_ancestry EAS EUR
    do
    
      mkdir -p ${RA}/${ancestry}/leads
      mkdir -p ${RA}/${ancestry}/summary_stats
    done
done

######################################### Summary stats

# Download data
sbatch --wrap="wget https://data.cyverse.org/dav-anon/iplant/home/kazuyoshiishigaki/ra_gwas/ra_gwas-10-28-2021.tar"
# Untar
tar -xvf ra_gwas-10-28-2021.tar

# Move summary stats into folders
mv 10-28-2021/Trans_all_auto-10-2021.txt.gz AllRA/multi_ancestry/summary_stats
mv 10-28-2021/EUR_all_auto-10-2021.txt.gz AllRA/EUR/summary_stats
mv 10-28-2021/EAS_all_auto-10-2021.txt.gz AllRA/EAS/summary_stats

mv 10-28-2021/Trans_seroposi_auto-10-2021.txt.gz SeropositiveRA/multi_ancestry/summary_stats
mv 10-28-2021/EUR_seroposi_auto-10-2021.txt.gz SeropositiveRA/EUR/summary_stats
mv 10-28-2021/EAS_seroposi_auto-10-2021.txt.gz SeropositiveRA/EAS/summary_stats

for RA in AllRA SeropositiveRA
do
  for ancestry in multi_ancestry EAS EUR
  do

    ## Reformat snp ID/ parse out chrom and position for lift over
    sbatch --mem=32G -t 72:00:00 --wrap="Rscript reformat_summary_stats.R ${RA} ${ancestry}"


    ## Liftover to hg38
    awk -F ',' '{print $1 "\t" $2 "\t" $2 "\t" $5}' ${RA}/${ancestry}/summary_stats/${RA}_${ancestry}.csv > ${RA}/${ancestry}/summary_stats/${RA}_${ancestry}_temp.bed

    sbatch -t 72:00:00 --mem=32G --wrap="Rscript remove_header.R ${RA} ${ancestry}"

    awk '{$3+=1}1' ${RA}/${ancestry}/summary_stats/${RA}_${ancestry}_temp_noheader.bed > ${RA}/${ancestry}/summary_stats/${RA}_${ancestry}_temp_end.bed

    sbatch -J ${RA}_${ancestry}_liftover --wrap="liftOver -bedPlus=1 ${RA}/${ancestry}/summary_stats/${RA}_${ancestry}_temp_end.bed \
      /proj/phanstiel_lab/References/genomes/hg19/hg19ToHg38.over.chain.gz ${RA}/${ancestry}/summary_stats/${RA}_${ancestry}_liftOver.bed \
       ${RA}/${ancestry}/summary_stats/${RA}_${ancestry}_liftOver_unmapped.bed"


    # Join liftover
    sbatch --mem=32G -t 6:00:00 -J ${RA}_${ancestry}_summaryLiftOver --wrap="Rscript join_liftOver_summary_stats.R ${RA} ${ancestry}"
  done
done





######################################### Leads

# 1) Grab variant allele frequencies from chrom vcfs for each ancestry
for chr in {1..22}
do
  # multi ancestry
  #sbatch -J multi_chr${chr} -t 72:00:00 --mem=8G --wrap="bcftools query -f '%ID %INFO/AF\n' /proj/phanstiel_lab/References/genomes/1000G/GRCh37/ALL/all_snps/1000G.GRCh38.20190312_chr${chr}_id.vcf.gz > /work/users/n/e/nekramer/RA_GWAS/AFs/chr${chr}_1000G_GRCh38_ALL_AF.txt"
  sbatch -J multi_chr${chr} -t 72:00:00 --mem=4G --wrap="plink --bfile /proj/phanstiel_lab/References/genomes/1000G/GRCh37/ALL/all_snps/1000G_phase3_chr${chr} --freq --out /work/users/n/e/nekramer/RA_GWAS/AFs/chr${chr}_1000G_GRCh37_ALL_AF"

  # EUR
  sbatch -J EUR_chr${chr} -t 72:00:00 --mem=8G --wrap="bcftools query -f '%ID %INFO/EUR_AF\n' /proj/phanstiel_lab/References/genomes/1000G/GRCh37/EUR/all_snps/EUR_1000G.GRCh37_chr${chr}.vcf.gz > /work/users/n/e/nekramer/RA_GWAS/AFs/chr${chr}_1000G_GRCh37_EUR_AF.txt"
  
  # EAS
  sbatch -J EAS_chr${chr} -t 72:00:00 --mem=8G --wrap="bcftools query -f '%ID %INFO/EAS_AF\n' /proj/phanstiel_lab/References/genomes/1000G/GRCh37/EAS/all_snps/EAS_1000G.GRCh37_chr${chr}_id.vcf.gz > /work/users/n/e/nekramer/RA_GWAS/AFs/chr${chr}_1000G_GRCh37_EAS_AF.txt"
done

# 2) Parse leads from supp table, separating top signals into ancestry/subtype called

sbatch -t 4:00:00 -J parseLeads --wrap="Rscript parseLeads.R Ishigaki_gwas_sup4.csv"

for RA in AllRA SeropositiveRA
do
    for ancestry in multi_ancestry EAS EUR
    do
      
      # 3) For each lead SNP, get LD buddies from corresponding 1000G ancestry, matching on hg19 variant IDs
      for snp in `awk -F ',' '{print $2}' ${RA}/${ancestry}/leads/${RA}_${ancestry}_leads.csv | awk 'NR!=1 {{print}}'`
      do
        # Get chrom
        chr=$(echo $snp | awk -F ':' '{print $1}')
        
        if [ "${ancestry}" == "multi_ancestry" ]; then
          sbatch -t 12:00:00 --mem=16G --wrap="plink --bfile /proj/phanstiel_lab/References/genomes/1000G/GRCh37/ALL/all_snps/1000G_phase3_chr${chr} --ld-snp ${snp} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out ${RA}/multi_ancestry/leads/${snp}_ld"
        elif [ "${ancestry}" == "EUR" ]; then
          sbatch -t 12:00:00 --mem=16G --wrap="plink --bfile /proj/phanstiel_lab/References/genomes/1000G/GRCh37/EUR/all_snps/EUR_1000G.GRCh37_chr${chr} --ld-snp ${snp} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out ${RA}/EUR/leads/${snp}_ld"
        else
          sbatch -t 12:00:00 --mem=16G --wrap="plink --bfile /proj/phanstiel_lab/References/genomes/1000G/GRCh37/EAS/all_snps/EAS_1000G.GRCh37_chr${chr}_id --ld-snp ${snp} --ld-window 200000 --ld-window-kb 1000 --ld-window-r2 0 --r2 --out ${RA}/EAS/leads/${snp}_ld"
        fi
        
        # Join LD buddies with lead data 
        sbatch -t 12:00:00 --mem=4G --wrap="Rscript reformat_LDbuddies.R ${snp} ${RA} ${ancestry} ${RA}/${ancestry}/leads/${snp}_ld.ld"
        
      done
      
      # 4) Concatenate all leads_ldbuddies, making sure to keep leads that weren't found in 1000G
      sbatch --wrap="Rscript concatenate_lead_lds.R ${RA} ${ancestry}"
      
      # 5) Lift over leads and LD buddies to hg38
      ## 5a) leads
      # 5ai) Create temporary BED-formatted file from reformatted leads/buddies
      awk -F ',' '{print $3 "\t" $4 "\t" $4 "\t" $1}' ${RA}/${ancestry}/leads/${RA}_${ancestry}_leads_ld.csv > ${RA}/${ancestry}/leads/${RA}_${ancestry}_leads_temp.bed
      # 5aii) Remove header
      sed -i '1d' ${RA}/${ancestry}/leads/${RA}_${ancestry}_leads_temp.bed
      
      # 5aiii) Add 1 bp for end position
      awk '{$3+=1}1' ${RA}/${ancestry}/leads/${RA}_${ancestry}_leads_temp.bed > ${RA}/${ancestry}/leads/${RA}_${ancestry}_leads_temp_end.bed
      
      # 5aiv) Lift over to hg38
      sbatch -J ${RA}_${ancestry}_liftover --wrap="liftOver -bedPlus=1 ${RA}/${ancestry}/leads/${RA}_${ancestry}_leads_temp_end.bed \
        /proj/phanstiel_lab/References/genomes/hg19/hg19ToHg38.over.chain.gz ${RA}/${ancestry}/leads/${RA}_${ancestry}_lead_liftOver.bed \
         ${RA}/${ancestry}/leads/${RA}_${ancestry}_leads_liftOver_unmapped.bed"

      ## 5b) buddies
      # 5bi) Create temporary BED-formatted file from reformatted leads/buddies
      awk -F ',' '{print $3 "\t" $11 "\t" $11 "\t" $12}' ${RA}/${ancestry}/leads/${RA}_${ancestry}_leads_ld.csv > ${RA}/${ancestry}/leads/${RA}_${ancestry}_buddies_temp.bed
      # 5bii) Remove header
      sed -i '1d' ${RA}/${ancestry}/leads/${RA}_${ancestry}_buddies_temp.bed

      # 5biii) Add 1 bp for end position
      awk '{$3+=1}1' ${RA}/${ancestry}/leads/${RA}_${ancestry}_buddies_temp.bed > ${RA}/${ancestry}/leads/${RA}_${ancestry}_buddies_temp_end.bed

      # 5biv) Lift over to hg38
      sbatch -J ${RA}_${ancestry}_liftover --wrap="liftOver -bedPlus=1 ${RA}/${ancestry}/leads/${RA}_${ancestry}_buddies_temp_end.bed \
        /proj/phanstiel_lab/References/genomes/hg19/hg19ToHg38.over.chain.gz ${RA}/${ancestry}/leads/${RA}_${ancestry}_buddies_liftOver.bed \
         ${RA}/${ancestry}/leads/${RA}_${ancestry}_buddies_liftOver_unmapped.bed"


    ## 5c) Join lift over leads and buddies
    sbatch -t 6:00:00 -J ${RA}_${ancestry}_joinLiftOver --wrap="Rscript join_liftOver.R ${RA} ${ancestry}"


    ## 6) Join summary stat info and allele freq info for buddies
    sbatch -t 72:00:00 --mem=32G -o ${RA}_${ancestry}_sumaf.out --wrap="Rscript join_af_sumstats.R ${RA} ${ancestry}"

    done
done







  


