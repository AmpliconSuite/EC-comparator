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

	DICT_COLS = {"#chr":"chromosome",
				 "bin":"chromosome",
				 "start":"start",
				 "end":"end",
				 "stop":"end",
				 "orientation":"stranded",
				 "strand":"stranded"}
	COLS_TOKEEP_SORTED = ["#chr_frag", "start_frag", "stop_frag", "circ_id", "cn", "strand"]
	DICT_COLS_TOKEEP = {"#chr":"#chr_frag",
						"chromosome": "#chr_frag",
						"start": "start_frag",
						"stop": "stop_frag",
						"end": "stop_frag",
						"circ_id": "circ_id",
						"cycle_id": "circ_id",
						"weight": "cn",
						"score": "cn",
						"proportions": "cn",
						"strand":"strand",
						"orientation":"strand"
						}

class DDT:
	EUCLIDIAN = 'euclidian'
	EUCLIDIAN_THRESHOLD = 3000

	EUCLIDIAN_NORM_L1 = 'euclidian_norm_l1'
	EUCLIDIAN_THRESHOLD_L1 = 0.3

	EUCLIDIAN_NORM_L2 = 'euclidian_norm_l2'
	EUCLIDIAN_THRESHOLD_L2 = 0.3

	AUC_TRIANGLE = 'auc_triangle'
	AUC_TRIANGLE_THRESHOLD = 500000

	AUC_TRAPEZE = 'auc_trapeze'
	AUC_TRAPEZE_THRESHOLD = 500000

	MATCH_ANGLE = 'match_angle'
	MATCH_ANGLE_THRESHOLD = 1.2

	RELATIVE_METRIC = 'match_score'
	RELATIVE_METRIC_THRESHOLD = 0.3

	UNMATCHED = 'unmatched'
	UNMATCHED_THRESHOLD = 10000

	HAMMING = 'hamming'
	HAMMING_NORM = 'hamming_norm'
	OVERLAP = 'overlap'

	COSINE_SIMILARITY = 'cosine_similarity'