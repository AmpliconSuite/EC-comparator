"""
Constants and other settings
"""

class HEADER:
	CIRC_ID = "circ_id"
	CHR = "#chr"
	START = "start"
	END = "end"
	STRAND = "strand"
	CN = "estimated_cn"
	CHR1 = "chr1"
	CHR2 = "chr2"
class PROPS:

	# visualization
	ALPHA = 'alpha'
	COLOR = 'color'
	MATCHED = {ALPHA:0.8}
	UNMATCHED = {ALPHA: 0.2,
			   COLOR: 'gray'}
	MAX_COST = 100000

	# compute similarity
	INF = 3000000000
