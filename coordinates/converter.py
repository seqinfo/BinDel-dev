bin_width = 100000

with open("locations.info.tsv", encoding = "UTF-8", mode = "r") as f, open("all_regions.bed",  encoding = "UTF-8", mode = "w") as out:
    header = f.readline()
    out.write(header.strip() + "\n")
    for line in f:
        line = line.strip().split("\t")
        chromosome, start, end, focus = line[0], int(line[1]), int(line[2]), line[3]        
    
        while start + bin_width < end:
            out.write(f"{chromosome}\t{start}\t{(min(start + bin_width - 1, end))}\t{focus}\n")
            start = start + bin_width
