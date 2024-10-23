"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 11:57 AM 08/10/23

Configs
"""

from AmpliconComparison.utils.utils import HEADER as h
from AmpliconComparison.utils.utils import DDT as d
from AmpliconComparison.metrics import compare as m


class Configs:
	DICT_DIST = {
		h.BREAKPOINT_DISTANCE: {
			# breakpoints matching distances
			d.EUCLIDIAN: {
				h.WEIGHT: 1,
				h.THRESHOLD: d.EUCLIDIAN_THRESHOLD,
				h.DEFINITION: m.euclidian_distance,
				h.ENABLE: True,
				d.JACCARD_DISTANCE: {
					h.WEIGHT: 1,
					h.ENABLE: True
				}
			},
			d.RELATIVE_METRIC: {
				h.WEIGHT: 1,
				h.THRESHOLD: d.RELATIVE_METRIC_THRESHOLD,
				h.DEFINITION: m.match_score,
				h.ENABLE: False,
				d.JACCARD_DISTANCE: {
					h.WEIGHT: 1,
					h.ENABLE: True
				}
			},
			d.JACCARD_DISTANCE: {
				h.WEIGHT: 1,
				h.DEFINITION: m.compute_breakpoint_distance,
				h.ENABLE: False
			}

		},
		h.GENOMIC_FOOTPRINT: {
			# genomic footprint distances
			d.HAMMING_NORM: {
				h.WEIGHT: 1,
				h.DEFINITION: m.get_hamming_score_norm,
				h.ENABLE: True
			},
			d.COSINE_DISTANCE: {
				h.WEIGHT: 1,
				h.DEFINITION: m.get_cosine_distance_cn,
				h.ENABLE: True
			},
			d.COPYNUMBER_JC: {
				h.WEIGHT: 1,
				h.DEFINITION: m.get_jc_distance_cn,
				h.ENABLE: True
			},
			d.FRAGMENTS_DISTANCE: {
				h.WEIGHT: 0.5,
				h.DEFINITION: m.get_overlap_fragments_weighted,
				h.ENABLE: True
			},
			d.CYCLES_DISTANCE: {
				h.WEIGHT: 0.5,
				h.DEFINITION: m.get_overlap_cycles_weighted,
				h.ENABLE: True
			}
		},
		h.OTHER_PARAMS: {}
	}
