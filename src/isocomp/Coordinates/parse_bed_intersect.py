import pandas as pd

def parse_bed_intersect(bed_path: str) -> pd.DataFrame
	"""
	Take as input the result of a bedtools multiple intersect. 
	return grouped datframe
	"""
	
	ranges = pd.read_csv(bed_path, sep="\t", header=None)
	
	grouped_ranges = ranges.groupby([0,1,2])

	return grouped_ranges
