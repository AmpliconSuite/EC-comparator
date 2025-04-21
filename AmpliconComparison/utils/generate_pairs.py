"""
Generate all pairs to compare
"""
import ast
import os
import argparse
import pandas as pd
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

def get_directory(i,df):
	return os.path.dirname(df.iloc[i,22])

def get_file(i,df):
	fname = os.path.basename(df.iloc[i,22])
	sample = fname.split("_AA_results")[0]
	amplicon = "amplicon" + str(int(df.iloc[i,1]))
	return f"{sample}_AA_results/{sample}_{amplicon}_cycles.bed", f"{sample}_{amplicon}"

def get_path(i,df, prefix_root):
	dirname1 = get_directory(i,df)
	filename1, sample1 = get_file(i,df)
	if prefix_root:
		return os.path.join(prefix_root,dirname1,filename1), sample1
	else:
		return os.path.join(dirname1,filename1), sample1

def overlap(i,j,df,genelist):
	genes1 = ast.literal_eval(df.iloc[i,5].replace("'","")) + ast.literal_eval(df.iloc[i,6].replace("'",""))
	genes2 = ast.literal_eval(df.iloc[j,5].replace("'","")) + ast.literal_eval(df.iloc[j,6].replace("'",""))
	genes_overlap = set(genes1) & set(genes2)
	
	# check if the two sample have at least one gene overlapping
	if not genes_overlap:
		return False, None, None
	
	# check if the samples contain the selected genes
	if genelist and not set(genes_overlap) & set(genelist):
		return False, None, None

	return	True, ",".join(genes1), ",".join(genes2)

def run_subprocess(i,df,project, prefix_root, root, genelist):
	
	df_out = pd.DataFrame(columns=["project",
								   "sample1","sample2",
								   "amplicon1_name","amplicon2_name",
								   "amplicon1_complexity","amplicon2_complexity",
								   "tissueorigin1","tissueorigin2",
								   "sampletype1","sampletype2",
								   "genes1","genes2",
								   "amplicon1_path","amplicon2_path","outdir"])
	
	for j in range(i,df.shape[0]):
		ok, genes1,genes2 = overlap(i,j,df,genelist)

		amplicon_path1, amplicon_name1 = get_path(i,df,prefix_root)
		amplicon_path2, amplicon_name2 = get_path(j,df,prefix_root)

		if ok:

			new_row = {
				"project": project,
				"reference": df.iloc[i,12],
				"sample1": df.iloc[i,0],
				"sample2": df.iloc[j,0],
				"amplicon1_classification": df.iloc[i,3],
				"amplicon2_classification": df.iloc[j,3],
				"amplicon1_complexity": df.iloc[j,7],
				"amplicon2_complexity": df.iloc[j,7],
				"tissueorigin1": df.iloc[i,13],
				"tissueorigin2": df.iloc[j,13],
				"sampletype1": df.iloc[i,14],
				"sampletype2": df.iloc[j,14],
				"genes1": genes1,
				"genes2": genes2,
				"amplicon1_path": amplicon_path1,
				"amplicon2_path": amplicon_path2,
				"outdir": os.path.join(root,f"{amplicon_name1}_{amplicon_name2}") if root else f"{amplicon_name1}_{amplicon_name2}"
			}
			df_out.loc[len(df_out)] = new_row
	return df_out

  
def generate_pairs(input,output,project,prefix_root=None,root=None,genelist=None):
	
	results = []
	df = pd.read_csv(input, header=0, sep=",",index_col=0)
	# keep only entries having an amplicon
	df = df[~df["AA amplicon number"].isna()]
 
	# Use ProcessPoolExecutor to run the subprocesses in parallel
	with ProcessPoolExecutor(max_workers=10) as executor:
		# Submit all commands to be run in parallel
		futures = {
			executor.submit(run_subprocess, i, df, project, prefix_root, root, genelist): i for i in range(0,df.shape[0])
		}

		# Process results as they complete
		for future in as_completed(futures):
			command = futures[future]
			try:
				result_df = future.result()
				results.append(result_df)
				print(f"Output from command '{command}'",result_df.shape)
			except Exception as exc:
				print(f"Command '{command}' generated an exception: {exc}")
	
	# Merge all resulting dataframes
	merged_df = pd.concat(results, ignore_index=True)
	return merged_df
	
	
	df_out.to_csv(output, sep="\t",header=True,index=None)

def parse_args():
	parser = argparse.ArgumentParser(
		description="Infer gene regulatory networks with significance estimation."
	)

	parser.add_argument(
		"--input",
		type=str,
		required=True,
		help="Path to input aggregated results (e.g., CSV or TSV)."
	)

	parser.add_argument(
		"--output", 
		type=str,
		required=True,
		help="Path to generated configuration file."
	)
	
	parser.add_argument(
		"--project",
		type=str,
		required=True,
		help="Project name."
	)
	
	parser.add_argument(
		"--prefix-root", 
		type=str,
		required=False,
		default=None,
		help="Root path to where the inputs are located."
	)
 
	parser.add_argument(
		"--root", 
		type=str,
		required=False,
		default=None,
		help="Root path to where the comparisons will be located."
	)

	parser.add_argument(
		"--filtered-genes",
		type=str,
		nargs="+",
		required=False,
		default=[],
		help="List of gene names to filter the results on. Provide as space-separated list or with repeated flag."
	)

	return parser.parse_args()


if __name__ == "__main__":
	args = parse_args()
	print("Input file:", args.input)
	print("Output file:", args.output)
	print("Prefix root path:", args.prefix_root)
	print("Root path:", args.root)
	print("Filtered genes:", args.filtered_genes)
	generate_pairs(args.input,
                args.output,
                args.project,
                args.prefix_root,
                args.root,
                args.filtered_genes)
	
