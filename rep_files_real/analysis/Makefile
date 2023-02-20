SCRIPTDIR=../scripts
READSFOLDER=../data/raw/
CLEANFOLDER=./cleaned_reads/

#script locations
NGMALIGNER=$(SCRIPTDIR)/ngm_aligner.sh
CREATECONSENSUS=$(SCRIPTDIR)/create_consensus.sh
CLEANREADS=$(SCRIPTDIR)/clean_reads.sh
DEGENERATESITES=$(SCRIPTDIR)/degenerate.py
TMPDIR=/scale_wlg_nobackup/filesets/nobackup/uoa03626/

#data
READS=$(shell find ${READSFOLDER} -name "*.fastq" -and -name "*R1*")
READNAMES=$(notdir $(READS))
CORREADS=$(subst .fastq,.cor.fq,$(READNAMES))
SLICEDNAMES=$(subst .cor.fq,.cor_sliced.fastq,$(CORREADS))
SLICEDREADS=$(addprefix $(CLEANFOLDER),$(SLICEDNAMES))

#Use to skip certain steps
#SLICEDREADS=$(shell find ${CLEANFOLDER}/sliced/ -name "*sliced.fastq" -or -name "*singletons.fastq") 

REF=../data/e_grandis/ref.fa

#######
#clean up reads and use them to iteratively clean up the reference
#######

#Steps 1-3 
$(SLICEDREADS) : $(READS)
	${CLEANREADS} -i ${READSFOLDER}

e_mel_3/e_mel_3.fa: e_mel_3/e_mel_3.bam e_mel_2/e_mel_2.fa
	${CREATECONSENSUS} -o $@ -r $(word 2, $^) -i $<

e_mel_2/e_mel_2.fa: e_mel_2/e_mel_2.bam e_mel_1/e_mel_1.fa
	${CREATECONSENSUS} -o $@ -r $(word 2, $^) -i $<

e_mel_1/e_mel_1.fa: e_mel_1/e_mel_1.bam ${REF}
	${CREATECONSENSUS} -o $@ -r $(word 2, $^) -i $<

e_mel_3/e_mel_3.bam: e_mel_2/e_mel_2.fa $(SLICEDREADS)
	 mkdir -p e_mel_3
	${NGMALIGNER} -s 0.3 -r $< -o $@ -i ${CLEANFOLDER}

e_mel_2/e_mel_2.bam: e_mel_1/e_mel_1.fa $(SLICEDREADS)
    
	mkdir -p e_mel_2
	${NGMALIGNER} -s 0.3 -r $< -o $@ -i ${CLEANFOLDER}

e_mel_1/e_mel_1.bam: ${REF} $(SLICEDREADS)
   $(info Executing Step Four Iteration One)
	mkdir -p e_mel_1
	${NGMALIGNER} -s 0.3 -r $< -o $@ -i ${CLEANFOLDER}


############
#make final alignment, recalibrate and dedup it according to GATK best practices workflow
############

file = alignment.bam alignment.dedup.bam e_mel_3/e_mel_3.dict e_mel_3/e_mel_3.fa alingment.firstcalls.vcf alingment.recal.table alignment.recal.bam
all: ${file}

#Step 6
alignment.bam: e_mel_3/e_mel_3.fa
	${NGMALIGNER} -r $< -o $@ -i ${READSFOLDER}

%.dedup.bam: %.bam
	picard -Xmx500g MarkDuplicates INPUT=$< OUTPUT=$@ METRICS_FILE=$@.metrics.txt CREATE_INDEX=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=$TMPDIR 

%.dict: %.fa
	picard -Xmx500g CreateSequenceDictionary REFERENCE=$< OUTPUT=$@ 

%.firstcalls.vcf: %.dedup.bam e_mel_3/e_mel_3.fa e_mel_3/e_mel_3.dict
	mkdir -p bqsr; parallel --halt 2 gatk HaplotypeCaller --heterozygosity 0.025 -R $(word 2, $^) -I $< -stand-call-conf 50 -ploidy 2 -L NC_052612.1 -L NC_052613.1 -L NC_052614.1 -L NC_052615.1 -L NC_052616.1 -L NC_052617.1 -L NC_052618.1 -L NC_052619.1 -L NC_052620.1 -L NC_052621.1 -L NC_052622.1 -O bqsr/{}_$@ ::: $(shell seq 1 11); picard SortVcf $(foreach num, $(shell seq 1 11), I=bqsr/$(num)_$@) O=$@ SEQUENCE_DICTIONARY=$(word 3, $^)

%.recal.table: %.firstcalls.vcf %.dedup.bam e_mel_3/e_mel_3.fa
	gatk BaseRecalibrator -I $(word 2, $^) -R $(word 3, $^) --known-sites $< -O $@

%.recal.bam: %.recal.table %.dedup.bam e_mel_3/e_mel_3.fa
	gatk ApplyBQSR -I $(word 2, $^) -R $(word 3, $^) --bqsr-recal-file $< -O $@

###
#lift over annotations from egrandis to new alignment
###

liftover/e_mel_1_repeatmask.gff3: liftover/egrandis_201_repeatmask.gff3 e_mel_1/e_mel_1.fa
	liftOver -gff $< e_mel_1.fa.chain $@ /dev/null

liftover/e_mel_2_repeatmask.gff3: liftover/e_mel_1_repeatmask.gff3 e_mel_2/e_mel_2.fa
	liftOver -gff $< e_mel_2.fa.chain $@ /dev/null

liftover/e_mel_3_repeatmask.gff3: liftover/e_mel_2_repeatmask.gff3 e_mel_3/e_mel_3.fa
	liftOver -gff $< e_mel_3.fa.chain $@ /dev/null

liftover/e_mel_3_repeatmask.bed: liftover/e_mel_3_repeatmask.gff3
	gff2bed <$< >$@

liftover/e_mel_1_genes.gff3: liftover/egrandis_201_genes.gff3 e_mel_1/e_mel_1.fa
	liftOver -gff $< e_mel_1.fa.chain $@ /dev/null

liftover/e_mel_2_genes.gff3: liftover/e_mel_1_genes.gff3 e_mel_2/e_mel_2.fa
	liftOver -gff $< e_mel_2.fa.chain $@ /dev/null

liftover/e_mel_3_genes.gff3: liftover/e_mel_2_genes.gff3 e_mel_3/e_mel_3.fa
	liftOver -gff $< e_mel_3.fa.chain $@ /dev/null

liftover/e_mel_3_genes.bed: liftover/e_mel_3_genes.gff3
	gff2bed <$< >$@

###
#get degenerate sites for estimating parameters
###

e_mel_3_cds.bed: liftover/e_mel_3_genes.gff3
	cat $< | awk '$$3 == "CDS" {print}' | gff2bed > $@

degenerate/e_mel_3_degenerate.bed: e_mel_3_cds.bed e_mel_3/e_mel_3.fa
	python2 ${DEGENERATESITES} -i $< -r $(word 2, $^) -o $@

# Reed's Rules
e_mel_3/e_mel_3.bed: e_mel_3/e_mel_3.fa.fai
	awk 'BEGIN {FS="\t"}; {print $$1 FS "0" FS $$2}' $< > $@

e_mel_3/e_mel_3_norepeats.bed: e_mel_3/e_mel_3.bed
	bedtools subtract -a $< -b liftover/e_mel_3_repeatmask.bed \
	| bedtools merge -d 10 \
	| awk '$$3-$$2 >= 200' \
	> $@

e_mel_3/e_mel_3_norepeats.regions: e_mel_3/e_mel_3_norepeats.bed
	awk '{print $$1 ":" $$2+1 "-" $$3'} < $< > $@

alignment_norepeats.ad:
	dng pileup -Q 13 -q 3 --max-dp 500 --min-dp 16 -o $@ -r @e_mel_3_norepeats.regions -f e_mel_3/e_mel_3.fa alignment.bam
