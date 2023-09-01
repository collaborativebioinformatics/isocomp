## Iso-Seq analysis
 
Demultiplexed hifi reads (Q20, single-molecule resolution) from lima [https://github.com/pacificbiosciences/barcoding/ https://github.com/pacificbiosciences/barcoding/] were 
further processed using Isoseq3 (v3.2.2) [https://github.com/PacificBiosciences/IsoSeq] refinement and clustering steps to generate 
the HQ (Full-length high quality i.e., predicted accuracy ≥Q20) transcripts [Table 1]. 
Next, we mapped these transcripts against GRCh38 (v33 p13) [2] using Minimap2 long read alignment tools [1] 
(v2.24-r1122; commands: minimap2 -t 8 -ax splice:hq -uf --secondary=no -C5 -O6,24 -B4 GRCh38.v33p13.primary_assembly.fa sample.polished.hq.fastq.gz). 
The table 2 shows the basic statistics of the alignment of each sample [HG002: NA27730, NA24385 and NA26105].
Next, we performed cDNA_cupcake v28.0.0 [https://github.com/Magdoll/cDNA_Cupcake] workflow to collapse the redundant isoforms from bam, 
followed by filtering the low counts isoforms by 10 and filter away 5' degraded isoforms as they are not biologically significant. Next, sqanti3 v5.0 [3] tool was
used to generate the final corrected fasta transcripts and gff [Table 3] along with the isoform classification reports. The external databases including reference data set of transcription start sites (refTSS), 
list of polyA motif, tappAS-annotation and Genecode GRCh38 annotation were utilized during the isoform classification by Sqanti3. Finally, IsoAnnotLite (v2.7.3) [https://isoannot.tappas.org/isoannot-lite/] 
analysis was performed to annotate the gtf file from sqanti3.

Table 1. Statistics of HQ fastq reads.
Samples	num_seqs	min_len	avg_len	max_len	N50	Q20(%)	Q30(%)
NA26105 (HG002.3)	205,590	51	3298.8	14887	3817	99.95	99.92
NA27730 (HG002.2)	205,884	52	3093.6	12524	3695	99.96	99.93
NA24385/HG002 (HG002.1)	411,349	50	2171.2	10767	2675	99.96	99.93

Table 2.  Statistics of isoforms mapped with Minimap2.
Samples	Total Reads	AlnPerc (%)	Mapped Reads	Unmapped reads
NA26105	MM2	99.04	203607	1983
NA27730	MM2	99.35	204553	1331
NA24385	MM2	99.65	409905	1444
 
Table 3. Basic Statistics of transcripts discovered by Sqanti3 filtered.
Sample	# of Gene	# of transcript	# of exon
NA26105	5433	7632	72940
NA27730	5459	7765	74906
NA24385/HG002	5656	9710	51196

Reference:
1.	Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18, September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191
2.	Adam Frankish and others, GENCODE 2021, Nucleic Acids Research, Volume 49, Issue D1, 8 January 2021, Pages D916–D923, https://doi.org/10.1093/nar/gkaa1087
3.	Pardo-Palacios FJ, Arzalluz-Luque A et al. SQANTI3: curation of long-read transcriptomes for accurate identification of known and novel isoforms. doi: 10.1101/2023.05.17.541248. PMID: 37398077; PMCID: PMC10312485.
![image](https://github.com/collaborativebioinformatics/isocomp/assets/41694905/60b9790a-3caf-499c-b799-bcd819ec0e6c)
