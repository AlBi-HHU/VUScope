with open(snakemake.output[0], "w") as outfile:
    for file in snakemake.input:
        with open(file, "r") as infile:
            outfile.write(infile.read() + "\n")