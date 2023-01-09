from os.path import join

configfile: "config/config.yaml"

# get list of accession numbers from accession number file
with open(config["accessionNumFile"], 'r') as f:
    accession_num_file = f.read()
    ACCESSION_NUMBERS = accession_num_file.split()

rule all:
    input:
# put which output files you want snakemake to produce here!
# for example, if you want to download the metagenomic FASTA files from NCBI's WGS db, put:
        expand(join(config["repositoryDir"],"{accession_num}"), accession_num=ACCESSION_NUMBERS),
        expand(join(config["dumpDir"],"{accession_num}.fa"),accession_num=ACCESSION_NUMBERS)

include:
# put which snakemake rule you want to run here!
# for example, if you want to run the snakemake rules to download the FASTA files above, put:
    "workflow/rules/download_metagenomic_data.smk"
