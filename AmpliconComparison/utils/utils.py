"""
Constants and other settings
"""
import json
import numpy as np

PACKAGE_NAME = 'AmpliconComparison'

class NpEncoder(json.JSONEncoder):
	def default(self, obj):
		if isinstance(obj, np.integer):
			return int(obj)
		if isinstance(obj, np.floating):
			return float(obj)
		if isinstance(obj, np.ndarray):
			return obj.tolist()
		return super(NpEncoder, self).default(obj)

class HEADER:
	CIRC_ID = "circ_id"
	CHR = "Chromosome"
	CHRH2 = "Chromosome2"
	START = "Start"
	END = "End"
	STRAND = "Stranded"
	CN = "estimated_cn"
	CHR1 = "chr1"
	CHR2 = "chr2"
	IDX1 = "idx1"
	IDX2 = "idx2"
	LEN = "len"
	S1 = "s1"
	S2 = "s2"
	TRACK  = "track"
	ISCYCLIC = "iscyclic"
	BINS = "bins"
	BIN_ENABLED = "bin_enabled"

	HEADER_SORTED = [CHR, START, END, CIRC_ID, CN, STRAND, ISCYCLIC]
	DICT_HEADER = {"#chr": CHR,
				   "chr": CHR,
				   "chromosome": CHR,
				   "start": START,
				   "stop": END,
				   "end": END,
				   "circ_id": CIRC_ID,
				   "cycle_id": CIRC_ID,
				   "weight": CN,
				   "score": CN,
				   "proportions": CN,
				   "estimated_cn": CN,
				   "cn": CN,
				   "strand": STRAND,
				   "orientation": STRAND}

	THRESHOLD = "threshold"
	WEIGHT = "weight"
	ENABLE = "enable"
	DEFINITION = "definition"

	BREAKPOINT_DISTANCE = 'breakpoint_distance'
	GENOMIC_FOOTPRINT = 'genomic_footprint'
	OTHER_PARAMS = 'other_params'
	CONFIGS = 'configs'
	DISTANCES = 'distances'
	DEFAULT = 'default'


class PROPS:
	# visualization
	ALPHA = 'alpha'
	COLOR = 'color'
	MATCHED = {ALPHA: 0.7}
	UNMATCHED = {ALPHA: 0.3,
				 COLOR: 'gray'}
	MAX_COST = 100000

	# compute similarity
	INF = 300000000

	DICT_COLS = {"#chr": "Chromosome",
				 "bin": "Chromosome",
				 "start": "Start",
				 "end": "End",
				 "stop": "End",
				 "orientation": "Stranded",
				 "strand": "Stranded"}
	COLS_TOKEEP_SORTED = ["#chr_frag", "start_frag", "stop_frag", "circ_id", "cn", "strand"]
	DICT_COLS_TOKEEP = {"#chr": "#chr_frag",
						"chromosome": "#chr_frag",
						"start": "start_frag",
						"stop": "stop_frag",
						"end": "stop_frag",
						"circ_id": "circ_id",
						"cycle_id": "circ_id",
						"weight": "cn",
						"score": "cn",
						"proportions": "cn",
						"strand": "strand",
						"orientation": "strand"
						}
	RNULL = "rnull"
	TNULL = "tnull"

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
	SIGMOID_THRESHOLD = 250

	HAMMING = 'cn_hamming_dist'
	HAMMING_NORM = 'cn_hamming_norm_dist'
	OVERLAP = 'overlap'
	OVERLAP_SCORE = 'overlapping_score'
	PROD = 'prod'
	LEN = 'len'
	TOTALLEN = 'totallen'

	COSINE_DISTANCE = 'cn_cos_dist'
	BREAKPOINT_SIMILARITY = 'breakpoint_similarity'
	BREAKPOINT_DISTANCE = 'breakpoint_dist'
	FRAGMENTS_DISTANCE = 'fragments_dist'
	CYCLES_DISTANCE = 'cycles_dist'
	JACCARD_DISTANCE = 'breakpoint_dist'
	COPYNUMBER_JC = 'cn_jc_dist'

	TOTAL_COST = 'total_cost'
	TOTAL_COST_DESCRIPTION = 'total_cost_description'

class OUTFILES:
	TOTAL_COST_PNG = 'total_cost.png'
	TOTAL_COST_TABLE = 'total_cost_table.png'
	COVERAGE_PROFILE_S1_TXT = 'coverage_profile_s1.txt'
	COVERAGE_PROFILE_S2_TXT = 'coverage_profile_s2.txt'
	COVERAGE_PROFILE_PNG = 'coverage_profile.png'
	COVERAGE_BREAKPOINTS_PROFILE = 'coverage_breakpoints_profile.png'
	BREAKPOINTS_PROFILE_S1_TXT = 'breakpoints_profile_s1.txt'
	BREAKPOINTS_PROFILE_S2_TXT = 'breakpoints_profile_s2.txt'
	METRICS_JSON = 'metrics.json'


