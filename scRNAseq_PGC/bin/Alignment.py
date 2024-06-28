'''
create genome index
STAR  --runMode genomeGenerate --runThreadN ... --genomeDir ./ --genomeFastaFiles /path/to/genome.fa  --sjdbGTFfile /path/to/genes.gtf

STAR --outSAMattributes All --outSAMtype BAM Unsorted --quantMode GeneCounts --readFilesCommand gunzip -c --runThreadN $NCPU
 --sjdbGTFfile $GTFFILE --outReadsUnmapped Fastx --outMultimapperOrder Random --genomeDir $GENOMEDIR
   --readFilesIn ${INPUTDIR}/${OUTPREFIX}_L001_R2_001.fastq.gz ${INPUTDIR}/${OUTPREFIX}_L001_R1_001.fastq.gz
    --outFileNamePrefix $OUTPREFIX --soloType CB_UMI_Simple --soloCBwhitelist $WHITELIST --soloUMIlen 12 --soloCBlen 16 --soloUMIstart 17

'''

