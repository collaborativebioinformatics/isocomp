# for each group in the grouped pandas dataframe, conduct the pairwise 
# comparison

# chr1	43969998	43973369	PB.547.14	0	+	GM27730	chr1	43969998	43973369	PB.580.12	0	+
# chr1	43969998	43973369	PB.547.14	0	+	GM27730	chr1	43969998	43973369	PB.580.17	0	+
# chr1	43969998	43973369	PB.547.14	0	+	HG002	chr1	43969998	43973369	PB.589.13	0	+
# chr1	43969998	43973369	PB.547.17	0	+	GM27730	chr1	43969998	43973369	PB.580.12	0	+
# chr1	43969998	43973369	PB.547.17	0	+	GM27730	chr1	43969998	43973369	PB.580.17	0	+
# chr1	43969998	43973369	PB.547.17	0	+	HG002	chr1	43969998	43973369	PB.589.13	0	+

# The 4th column (count starting at 1) is the "reference" isoform. Note that 
# there are two isoforms from the "reference" in this region, PB.547.14 and 
# PB.547.14.
# Column 7 describes the "query" sample and column 11 describes the "query"
# transcript.

# We must know the sample name of the "reference" sample (this is input)
# Using the "sample" + "isoform_name", we can extract the sequence from 
# the corresponding fasta files (input)
# We can do the same thing using the query "sample" and query "transcript"

# We can group the comparisons by reference isoform. A query sequence that 
# has no match against any of the reference query is an instance of a isoform 
# that has "equivalent" coordinates, but possesses a sequence variant