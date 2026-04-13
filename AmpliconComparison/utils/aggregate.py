import argparse
import pandas as pd
import json
import os
from pathlib import Path

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Filter a file and read JSON data.")
    parser.add_argument("--input-file", required=True, help="Path to input pairwise comparison specification.")
    parser.add_argument("--output-file", required=True, help="Path to file containing the aggregated data.")
    args = parser.parse_args()

    # Read the input file into a DataFrame (no header)
    try:
        df = pd.read_csv(args.input_file, sep="\t", header=0)
    except Exception as e:
        print(f"Error reading the input file: {e}")
        return

    # Open the output file for writing
    try:
        with open(args.output_file, "w") as g:
            
            g.write("\t".join(map(str,df.columns.tolist())) + "\t" + "total_cost\t" + "cn_hamming_norm_dist\t" + "cn_cos_dist\t" + "cn_jc_dist\t" + "breakpoint_dist\t" + "fragments_dist\t" + "cycles_dist\n")

            for index, value in enumerate(df["outdir"]):
                
                data_str = ""

                file_path=Path(os.path.join(value,"metrics.json"))

                if file_path.is_file():
                    # Read the JSON file
                    with open(os.path.join(value, "metrics.json"), 'r') as file:
                        data = json.load(file)
                
                    # add total cost
                    data_str += str(data["distances"]["total_cost"]) + "\t"
                
                    # add metrics contributing to the total cost
                    data_str += str(data["distances"]["cn_hamming_norm_dist"]) + "\t"
                    data_str += str(data["distances"]["cn_cos_dist"]) + "\t"
                    data_str += str(data["distances"]["cn_jc_dist"]) + "\t"
                    data_str += str(data["distances"]["breakpoint_dist"]) + "\t"
                    data_str += str(data["distances"]["fragments_dist"]) + "\t"
                    data_str += str(data["distances"]["cycles_dist"]) + "\t"

                    # Write the value to the output file
                    g.write("\t".join(map(str, df.iloc[index].values)) + "\t" + data_str + "\n")
                else:
                    print(f"Not found: {value}")  # Optional: Print to console
    except Exception as e:
        print(f"Error writing to the output file: {e}")

if __name__ == "__main__":
    main()