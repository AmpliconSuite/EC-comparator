"""
Constants and other settings
"""
import json
import numpy as np

PACKAGE_NAME = 'AmpliconComparison'

def get_weight_distance(key, d, dict_metrics):

	if key == HEADER.GENOMIC_FOOTPRINT:
		# check if distance was enabled
		if d in dict_metrics[HEADER.CONFIGS][key] and dict_metrics[HEADER.CONFIGS][key][d][HEADER.ENABLE]:
			# distance is enabled
			weight = dict_metrics[HEADER.CONFIGS][key][d][HEADER.WEIGHT]
			return weight, d

	elif key == HEADER.BREAKPOINT_DISTANCE:
		def_dist = dict_metrics[HEADER.CONFIGS][key][HEADER.DEFAULT]
		if def_dist in dict_metrics[HEADER.CONFIGS][key] and dict_metrics[HEADER.CONFIGS][key][def_dist][HEADER.ENABLE]:
			weight = dict_metrics[HEADER.CONFIGS][key][def_dist][HEADER.WEIGHT]
			return weight, def_dist

	return None, None


def get_value_distance(key, d, dict_metrics):

	if key == HEADER.GENOMIC_FOOTPRINT:
		if d in dict_metrics[HEADER.CONFIGS][key] and dict_metrics[HEADER.CONFIGS][key][d][HEADER.ENABLE]:
			val = dict_metrics[HEADER.DISTANCES][d]
			return val, d

	elif key == HEADER.BREAKPOINT_DISTANCE:
		def_dist = dict_metrics[HEADER.CONFIGS][key][HEADER.DEFAULT]
		if def_dist in dict_metrics[HEADER.CONFIGS][key] and dict_metrics[HEADER.CONFIGS][key][def_dist][HEADER.ENABLE]:
			val = dict_metrics[HEADER.DISTANCES][HEADER.BREAKPOINT_DISTANCE]
			return val, HEADER.BREAKPOINT_DISTANCE

	return None, None



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
	FRAG_ID = "Frag_id"
	BINS = "bins"
	BIN_ENABLED = "bin_enabled"

	HEADER_SORTED = [CHR, START, END, CIRC_ID, CN, STRAND, ISCYCLIC, FRAG_ID]
	DICT_HEADER = {"#chr": CHR,
				   "chr": CHR,
				   "chromosome": CHR,
				   "start": START,
				   "stop": END,
				   "end": END,
				   "circ_id": CIRC_ID,
				   "cycle_id": CIRC_ID,
				   "structure": CIRC_ID,
				   "fragment": FRAG_ID,
				   "fragment_id": FRAG_ID,
				   "weight": CN,
				   "score": CN,
				   "proportions": CN,
				   "estimated_cn": CN,
				   "estimated_proportions":CN,
				   "cn": CN,
				   "direction": STRAND,
				   "strand": STRAND,
				   "orientation": STRAND,
	   			   "iscyclic":ISCYCLIC,
			 	   "Iscyclic":ISCYCLIC}

	THRESHOLD = "threshold"
	WEIGHT = "weight"
	ENABLE = "enable"
	DEFINITION = "definition"

	BREAKPOINT_DISTANCE = 'breakpoint_dist'
	BREAKPOINT_DISTANCE_CALCULATION = 'breakpoint_dist_calc'
	GENOMIC_FOOTPRINT = 'genomic_footprint'
	OTHER_PARAMS = 'other_params'
	CONFIGS = 'configs'
	DISTANCES = 'distances'
	DEFAULT = 'default'


class PROPS:
	# visualization
	ALPHA = 'alpha'
	COLOR = 'color'
	MATCHED = {ALPHA: 0.9}
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
						"estimated_proportions":"cn",
						"strand": "strand",
						"orientation": "strand"
						}
	RNULL = "rnull"
	TNULL = "tnull"

class DDT:
	
	# breakpoints candidates selection using Gaussian distribution as cost function
	GAUSSIAN = 'gaussian'
	GAUSSIAN_DEF_SIGMA = 'sigma'
	GAUSSIAN_SIGMA = 500
	GAUSSIAN_DEF_AMPL = 'amplitude'
	GAUSSIAN_AMPL = 1
 	   
	# breakpoints candidates selection using Euclidian distance as cost function
	EUCLIDIAN = 'euclidian'
	EUCLIDIAN_THRESHOLD = 3000
	
	EUCLIDIAN_WEIGHTED = 'euclidian_weighted'
	EUCLIDIAN_WEIGHTED_THRESHOLD = 3000

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

	UNMATCHED_DISTANCE = 'unmatched_distance'
	UNMATCHED_THRESHOLD = 'umatched_treshold'
	MANHATTAN = 'manhattan'
	MANHATTAN_THRESHOLD = 1000
	SIGMOID_THRESHOLD = 1000

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
 
	COVERAGE = 'coverage'
	MEAN_COVERAGE = 'mean_coverage'
	BP_MATCH_UNWEIGHTED = 'breakpoint_match_unweighted'
	BP_MATCH_CN_WEIGHTED = 'breakpoint_match_cn_weighted'
	BP_MATCH_CN_WEIGHTED_AVG = 'breakpoint_match_cn_weighted_avg'
	# in combination with gaussian  breakpoint matching
	BP_MATCH_UNWEIGHTED_CONFIDENCE = 'breakpoint_match_unweighted_confidence'
	BP_MATCH_CN_WEIGHTED_CONFIDENCE = 'breakpoint_match_cn_weighted_confidence'
	BP_MATCH_CN_WEIGHTED_AVG_CONFIDENCE = 'breakpoint_match_cn_weighted_avg_confidence'
	BP_GAUSSIAN_CONFIDENCE_UNWEIGHTED = 'breakpoint_gaussian_confidence_unweighted'
	BP_GAUSSIAN_CONFIDENCE_CN_WEIGHTED = 'breakpoint_gaussian_confidence_cn_weighted'
	
	BP_OPTIONS = [BP_MATCH_UNWEIGHTED, BP_MATCH_CN_WEIGHTED, BP_MATCH_CN_WEIGHTED_AVG,
			   	BP_MATCH_UNWEIGHTED_CONFIDENCE,BP_MATCH_CN_WEIGHTED_CONFIDENCE,BP_MATCH_CN_WEIGHTED_AVG_CONFIDENCE,
				BP_GAUSSIAN_CONFIDENCE_UNWEIGHTED,BP_GAUSSIAN_CONFIDENCE_CN_WEIGHTED]
	BP_CANDIDATES_SELECTION = [MANHATTAN, GAUSSIAN, RELATIVE_METRIC]
	BP_COST_FUNCTION = [EUCLIDIAN, GAUSSIAN, RELATIVE_METRIC]
	
	MATCH_NONLINEAR = 'match_nonlinear'
	MATCH_NONLINEAR_F1 = 'neg_exp'
	MATCH_NONLINEAR_F2 = 'hssp_curve'
	BP_OPTIONS_NONLINEAR = [MATCH_NONLINEAR_F1,MATCH_NONLINEAR_F2]
	MATCH_NONLINEAR_F1_THRESHOLD = 1000
	MATCH_NONLINEAR_F2_THRESHOLD = 1000
 
	RENAME = {COSINE_DISTANCE:'d2',
			 COPYNUMBER_JC:'d3',
             HAMMING_NORM: 'd1',
             FRAGMENTS_DISTANCE: 'd4',
             CYCLES_DISTANCE: 'd5',
             BREAKPOINT_DISTANCE: 'd6',
             JACCARD_DISTANCE:'d6'
             }


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
	MATCHED_BREAKPOINTS_TXT = 'breakpoints_matched.txt'


