import os
from argparse import ArgumentParser

parser: ArgumentParser = ArgumentParser()

parser.add_argument("-i", "--input", dest="in_file", help="input file")
parser.add_argument("-o", "--output", dest="out", help="output file name")

args = parser.parse_args()

counts = "focus\tHMM\tlength\n"

with open(os.path.abspath(args.in_file), mode="r", encoding="UTF-8") as f:
    f.readline()
    prev_focus, prev_HMM, prev_start, current_state_length = None, None, None, 0

    for line in f:
        line = line.strip().split("\t")
        focus, HMM, _ = line[0], line[1], int(line[2])

        if prev_focus is None:
            prev_focus = focus
            prev_HMM = HMM

        if prev_focus == focus and prev_HMM == HMM:
            current_state_length += 1
        else:
            counts += f"{prev_focus}\t{prev_HMM}\t{current_state_length}\n"
            current_state_length = 0
            prev_focus = focus
            prev_HMM = HMM

with open(os.path.abspath(args.out), mode="w", encoding="UTF-8") as out:
    out.write(counts)
