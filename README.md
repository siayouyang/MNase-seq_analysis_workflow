# MNase-seq_analysis_workflow
MNase-seq workflow for "Structure of the ISW1a complex bound to the dinucleosome"(2023,NSMB)

1.	unzip：
nohup gunzip wt1_1.fq.gz &

2.	remove adapter:
nohup java -jar /data02/octopus/Octopus-toolkit/Octopus-toolkit/Tools/Trimmomatic/trimmomatic.jar PE wt1_1.fq wt1_2.fq -baseout wt1 ILLUMINACLIP:/data02/octopus/Octopus-toolkit/Octopus-toolkit/Tools/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:1:true SLIDINGWINDOW:4:15 &

3.	Bowtie:
nohup bowtie2 -p 5 --no-mixed --no-discordant --no-unal -x /data01/hpj/reference_genom/sacCer3/genome -1 wt1_1P -2 wt1_2P -S wt1.sam &

4.	Fragment length distribution:
nohup /data02/xyy/bash_scripts/sort.sh wt1.sam > wt1.txt &

5.	Select 120-180bp fragments, remove rRNA region fragments:
nohup python /data02/xyy/python_scripts/SAMFilter.py -i wt1.sam --minLength 120 --maxLength 180 -o wt1.filtered.sam &

6.	Remove PCR duplicates:
nohup /data03/xyy/bash_scripts/removeduplicate.sh wt1.filtered.sam &

7.	Convert BAM to SAM:
nohup samtools view -h wt1.filtered.samrd.bam > wt1.filtered.samrd.sam

8.	Count reads:
python /data03/xyy/python_scripts/CountRandSelect.py -i wt1.filtered.samrd.sam -c &

9.	Random select reads（downsampling）:
nohup python /data03/xyy/python_scripts/CountRandSelect.py -i wt1.filtered.samrd.sam -n 13480875 &

10.	Find dyad position（slow）:
nohup python /data03/xyy/python_scripts/SAMtoPosition.py --rawWig --smoothWig -p 20 -i wt1.filtered.samrd.sam.randsel.sam &

11.	Wig to bigWig（optional）:
nohup /data02/xyy/python_scripts/wigToBigWig wt1.filtered.samrd.sam.randsel.sam.smooth.wig /data02/xyy/python_scripts/sacCer3.chrom.sizes wt1.filtered.samrd.sam.randsel.sam.smooth.bw

12.	Compare position shift（execute one at a time）:
nohup python /data03/xyy/python_scripts/CompareNucleosome.py --input1 wt1.filtered.samrd.sam.randsel.sam.info.txt --input2 wt2.filtered.samrd.sam.randsel.sam.info.txt --output1 wt1_wt2 --output2 wt2_wt1 -p 20 -g /data04/xyy/python_scripts/clark_genome_2014_TSS_TTS.bed --significantShiftGenes &

13.	Get overlapping genes (.txt) from http://www.biovenn.nl/

14.	make_overlap_bed.py（manually change the input and output of .bed, .txt inside the scripts）

15.	use the overlappinggenes (.bed) as input to run step 12.（set --pvalue 1）:
nohup python /data03/xyy/python_scripts/CompareNucleosome.py --input1 wt1.filtered.samrd.sam.randsel.sam.info.txt --input2 wt2.filtered.samrd.sam.randsel.sam.info.txt --output1 wt1_wt2_1108genes --output2 wt2_wt1_1108genes -p 20 -g /data04/xyy/python_scripts/overlap_1108genes.bed --significantShiftGenes --pvalue 1 &

17.	composite plot(execute one at a time, fast):
nohup python /data04/xyy/python_scripts/SelectWindowCompositePlot.py -i wt1.filtered.samrd.sam.randsel.sam.smooth.wig -p 20 -o wt1_composite -g ../1108genes_sigshift_bed/wt1_wt2_1108genes.TSS4nuc.sigshift.bed -b 100 -a 600
