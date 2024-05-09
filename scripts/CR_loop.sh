
for S in mRNA_hiPSC mRNA_hiPSC_CM;
do
	/20tb/ratto/Bertero/ForPaper/5/cellranger-7.0.0/cellranger count --id $S --fastqs=/20tb/ratto/Bertero/ForPaper/5/201116,/20tb/ratto/Bertero/ForPaper/5/201117 --sample="$S" --transcriptome /20tb/ratto/Bertero/ForPaper/5/cellranger-7.0.0/refdata-gex-GRCh38-2020-A 
done

mRNA_D0
for S in mRNA_D2 mRNA_D6 mRNA_D12;
do
  echo $S
  /20tb/ratto/Bertero/ForPaper/5/cellranger-7.0.0/cellranger count --id $S --fastqs=/20tb/ratto/Bertero/ForPaper/7/fastq --sample="$S" --transcriptome /20tb/ratto/Bertero/ForPaper/5/cellranger-7.0.0/refdata-gex-GRCh38-2020-A 
done


/20tb/ratto/Bertero/ForPaper/5/cellranger-7.0.0/cellranger aggr --id=H001AS5e7 --csv=/20tb/ratto/Bertero/ForPaper/7/aggr_all_config.csv --normalize=none 

