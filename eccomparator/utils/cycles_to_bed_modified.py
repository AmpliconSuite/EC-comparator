#!/usr/bin/env python

import argparse
from collections import defaultdict
import os

from intervaltree import IntervalTree
import pandas as pd

# Write cycles files into merged beds

# read a cycles file into a dictionary of interval tree dictionaries
def read_cycles_file(fname):
    # all_cycles maps cycle_num -> chromosome -> IntervalTree
    seglookup = {}
    all_cycles_ivald = {}
    all_cycles = {}
    with open(fname) as infile:
        for line in infile:
            if line.startswith("Segment"):
                fields = line.rstrip().rsplit()
                segnum = int(fields[1])

                chr_ = fields[2]
                start_ = int(fields[3])
                end_ = int(fields[4])

                # avoid 0 size intervals
                end_ = end_ if start_ != end_ else end_+1
                dtup = (chr_, start_, end_)

                seglookup[segnum] = dtup

            elif line.startswith(("Cycle")):
                cycle_ivalt = defaultdict(IntervalTree)
                slist = []
                fields = line.rstrip().rsplit(';')
                cd = {x.rsplit("=")[0]: x.rsplit("=")[1] for x in fields}
                cnum = int(cd["Cycle"])
                segstring = cd["Segments"]
                copy_count = cd["Copy_count"]
                seglist = [(int(x[:-1]), x[-1]) for x in segstring.rsplit(",")]
                for x in seglist:
                    if x[0] == 0:
                        continue

                    dtup = seglookup[x[0]]
                    slist.append(dtup + (x[1], copy_count))
                    cycle_ivalt[dtup[0]].addi(dtup[1], dtup[2], (x[1], copy_count))  # duplicates will be overwritten

                all_cycles[cnum] = slist
                all_cycles_ivald[cnum] = cycle_ivalt

    return all_cycles_ivald, all_cycles

def write_bed(prefix, merged_intervals):
    with open(prefix + ".bed", 'w') as outfile:
        for i in merged_intervals:
            oline = "\t".join([str(x) for x in i]) + "\n"
            outfile.write(oline)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert cycles file to a series of bed files")
    parser.add_argument("-c", "--cycles", help="Cycles file", type=str, required=True)
    parser.add_argument("-o", "--output", help="Output file", type=str, required=True)
    args = parser.parse_args()

    print("In: ",args.cycles)
    print("Out: ",args.output)
    print()

    a_cycles_ivald, a_cycles_list = read_cycles_file(args.cycles)

    df_output = pd.DataFrame(columns=["#chr",	"start",	"end",	"strand",	"circ_id",	"estimated_cn"])
    for cnum, curr_cycle in a_cycles_ivald.items():
        for t in a_cycles_list[cnum]:
            data = {"#chr":t[0],
                    "start":t[1],
                    "end":t[2],
                    "strand":t[3],
                    "circ_id":str(cnum),
                    "estimated_cn":t[4]}
            df_output.loc[len(df_output)] = data

    df_output.to_csv(args.output, header=True, index=False, sep="\t")
