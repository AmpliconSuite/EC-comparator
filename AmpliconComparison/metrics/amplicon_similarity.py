#!/usr/bin/env python3
"""
This code is from https://github.com/AmpliconSuite/AmpliconClassifier/blob/main/ampclasslib/amplicon_similarity.py
It includes the Amplicon similarity score from Luebeck et al 2023, Nature 
https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-05937-5/MediaObjects/41586_2023_5937_MOESM1_ESM.pdf

But adapted to work with BED files only.

"""
import os
import argparse
from collections import defaultdict
import compare

add_chr_tag = False
bpweight = 0.75
d = 250
cn_cut = 4.5
min_de = 1
min_de_size = 2500


class Breakpoint(object):
	def __init__(self, lchrom, lpos, rchrom, rpos, cn):
		self.lchrom = lchrom
		self.lpos = lpos
		self.rchrom = rchrom
		self.rpos = rpos
		self.cn = cn

	def to_string(self):
		return self.lchrom + ":" + str(self.lpos) + " | " + self.rchrom + ":" + str(self.rpos) + "\t" + str(self.cn)

	def d_similar(self, bp2, d):
		bp2_chrom_set = {bp2.lchrom, bp2.rchrom}
		if self.lchrom not in bp2_chrom_set or self.rchrom not in bp2_chrom_set:
			return False

		sbp1 = sorted([(self.lchrom, self.lpos), (self.rchrom, self.rpos)])
		sbp2 = sorted([(bp2.lchrom, bp2.lpos), (bp2.rchrom, bp2.rpos)])

		if sbp1[0][0] == sbp2[0][0] and sbp1[1][0] == sbp2[1][0]:
			if abs(sbp1[0][1] - sbp2[0][1]) + abs(sbp1[1][1] - sbp2[1][1]) < d:
				return True

		return False


def amps_overlap(gdoi1, gdoi2):
	for c1, it1 in gdoi1.items():
		it2 = gdoi2[c1]
		if not it2 or not it1:
			continue

		for int1 in it1:
			if int1.data is None or int1.data > cn_cut:
				for int2 in it2[int1.begin:int1.end]:
					if int2.data is None or int2.data > cn_cut:
						return True

	return False


def nucSimilarity(g1, g2):
	g1_bp = 0
	obp = 0
	for c1, it1 in g1.items():
		for t1 in it1:
			g1_bp += (t1.end - t1.begin)
			for t2 in g2[c1][t1.begin:t1.end]:
				obp += (min(t1.end, t2.end) - max(t1.begin, t2.begin))

	try:
		nS = obp/g1_bp
		return nS

	except ZeroDivisionError:
		return 0


def bp_dist(bplist1, bplist2, d):
	intersectCount = 0.0
	for bp1 in bplist1:
		for bp2 in bplist2:
			if bp1.d_similar(bp2,d):
				intersectCount += 1
				break

	tE = max(len(bplist1), 1)
	return intersectCount/tE


def jaccard_sim_seq(g1, g2):
	g1_bp = 0
	g2_bp = 0
	obp = 0
	for c1, it1 in g1.items():
		for t1 in it1:
			g1_bp += (t1.end - t1.begin)
			for t2 in g2[c1][t1.begin:t1.end]:
				obp += (min(t1.end, t2.end) - max(t1.begin, t2.begin))

	for c2, it2 in g2.items():
		for t2 in it2:
			 g2_bp += (t2.end - t2.begin)

	jsim = obp/(g1_bp + g2_bp - obp)
	return jsim, obp, g1_bp, g2_bp


def compute_num_shared_bps(bplist1, bplist2, d):
	if not bplist1 and not bplist1:
		return 0

	intersectCount = 0.0
	for bp1 in bplist1:
		for bp2 in bplist2:
			if bp1.d_similar(bp2, d):
				intersectCount += 1
				break

	return intersectCount


def jaccard_sim_bp(bplist1, bplist2, d):
	if not bplist1 and not bplist1:
		return 0, 0

	intersectCount = compute_num_shared_bps(bplist1, bplist2, d)
	jsim = intersectCount/(len(bplist1) + len(bplist2) - intersectCount)
	return jsim, intersectCount


def asymmetric_score(bplist1, bplist2, st1, st2, d):
	ns = nucSimilarity(st1, st2)
	bpd = bp_dist(bplist1, bplist2, d)
	S = (1 - bpweight) * ns + (bpweight) * bpd
	return S, ns, bpd


def symScore(bplist_a, bplist_b, st_a, st_b, d):
	as1, nS1, bpd1 = asymmetric_score(bplist_a, bplist_b, st_a, st_b, d)
	as2, nS2, bpd2 = asymmetric_score(bplist_b, bplist_a, st_b, st_a, d)
	return (as1 + as2)/2, [as1, as2, nS1, nS2, bpd1, bpd2]

def parse(df1, df2, subset_ivald, cn_cut, min_de=1):
	"""
	Adapted to read the breakpoint graph from BED like format.
	"""
	bps = []
	no_subset = (len(subset_ivald) == 0)
	cn_cut_0 = (cn_cut == 0 and no_subset)
	keepAll = (no_subset or cn_cut_0)
	segTree = defaultdict(IntervalTree)
	
	bins, chrlist = compare.bin_genome(df1, df2, margin_size=0)

	# as pyranges objects
	df_t_pr, df_r_pr, df_bins_pr = compare.read_pyranges(df1, df2, bins)
	# as dataframe
	df_ = compare.overlap_vector(df_t_pr, df_r_pr, df_bins_pr)

	# 2. Compute hamming distance and others
	df_, h = compare.get_hamming_score(df_)

	# 3. Compute copy-number similarity
	cn_profile_t = get_feature_cn(df1, df_)
	cn_profile_r = get_feature_cn(df2, df_)
 
	
	
	
	with open(bpgf) as infile:
		
		for line in infile:
			# if line.startswith("discordant") or line.startswith("concordant"):
			if line.startswith("discordant"):
				fields = line.rstrip().rsplit()
				l, r = fields[1].rsplit("->")
				lchrom, lpos = l[:-1].rsplit(":")
				rchrom, rpos = r[:-1].rsplit(":")
				ldir, rdir = l[-1], r[-1]
				if add_chr_tag and not lchrom.startswith('chr'):
					lchrom = "chr" + lchrom
					rchrom = "chr" + rchrom

				lpos, rpos = int(lpos), int(rpos)
				if lcD[lchrom][lpos] or lcD[rchrom][rpos]:
					continue

				elif not keepAll and not (subset_ivald[lchrom].overlaps(lpos) and subset_ivald[rchrom].overlaps(rpos)):
					continue

				elif keepAll and not (segTree[lchrom].overlaps(lpos) or segTree[rchrom].overlaps(rpos)):
					continue

				elif lchrom == rchrom and abs(lpos - rpos) < min_de_size and ldir != rdir:
					continue

				cn = float(fields[2])
				currBP = breakpoint(lchrom, lpos, rchrom, rpos, cn)
				bps.append(currBP)

			elif line.startswith("sequence"):
				fields = line.rstrip().rsplit()
				lchrom, lpos = fields[1].rsplit(":")
				lpos = int(lpos[:-1])
				rchrom, rpos = fields[2].rsplit(":")
				rpos = int(rpos[:-1])+1

				if add_chr_tag and not lchrom.startswith('chr'):
					lchrom = 'chr' + lchrom
					rchrom = 'chr' + rchrom

				if lcD[lchrom][lpos:rpos] or cg5D[lchrom][lpos:rpos]:
					continue

				elif not keepAll and not subset_ivald[lchrom].overlaps(lpos, rpos):
					continue

				cn = float(fields[3])
				if cn > cn_cut:
					segTree[lchrom].addi(lpos, rpos, cn)

	if not bps and min_de > 0:
		return bps, defaultdict(IntervalTree)

	return bps, segTree

def compute_similarity(file1, file2, tempdir):

	df1, df2 = read_input(file1, file2, tempdir, min_cn=0)


	s, featList = symScore(bplist_a, bplist_b, st_a, st_b, d)
	jaccard_seq, amp_olap_len, amp_a_len, amp_b_len = jaccard_sim_seq(st_a, st_b)
	jaccard_bp, num_shared_bps = jaccard_sim_bp(bplist_a, bplist_b, d)
	featList.extend([jaccard_seq, jaccard_bp, num_shared_bps, len(bplist_a), len(bplist_b), amp_olap_len,
						 amp_a_len, amp_b_len])
	return [file1, file2, s] + featList
	
	
	return


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Compute similarity between overlapping AA amplicons")
	parser.add_argument("-a", help="Amplicon 1 /path/to/sample1_amplicon1_cycles.bed", required=True)
	parser.add_argument("-b", help="Amplicon 2 /path/to/sample2_amplicon1_cycles.bed", required=True)
	parser.add_argument("-o", help="Output dir")
	args = parser.parse_args()
	
	os.makedirs(args.o, exist_ok=True)
   
	with open(args.o + "/amplicon_similarity_scores.tsv", 'w') as outfile:
		# [as1, as2, nS1, nS2, bpd1, bpd2]
		outfile.write("Amp1\tAmp2\tSimilarityScore\tAsymmetricScore1\tAsymmetricScore2\t"
					  "GenomicSegmentScore1\tGenomicSegmentScore2\t"
					  "BreakpointScore1\tBreakpointScore2\t"
					  "JaccardGenomicSegment\tJaccardBreakpoint\t"
					  "NumSharedBPs\tAmp1NumBPs\tAmp2NumBPs\tAmpOverlapLen\tAmp1AmpLen\tAmp2AmpLen\n")
		outx = compute_similarity(as1, as2, bp1, bp2, args.o)


