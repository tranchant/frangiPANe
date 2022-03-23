#!/usr/bin/env python3

# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
# Intellectual property belongs to IRD, ... and SouthGreen development platform
# Written by Clothilde Chenal and Christine Tranchant-Dubreuil
# Copyright 2021

"""

    ....

    :Example:

    >>> from tools import frangiPANe

"""
import os
from os import path
import subprocess
import logging

import pandas as pd
import panel as pn

import shutil

from IPython.display import display, HTML
from tools.jupyter import display_alert


def make_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        at = "success"
        text = f"Directory {path} created"

    else:
        at = "warning"
        text = f"""### {at}
<hr>
Directory {path} already existed"""

    display_alert(text, at)


def read_group_file(file, logger):
    df_group = pd.read_csv(file, names=["sample", "Species"], index_col=False, sep="\t")
    id_dict = {}
    for line in open(file).readlines():
        field = line[:-1].split()  # [0] id, [1] group
        id_dict[field[0]] = field[1]

    logger.info(f"CHECKING GROUP FILE  OK")
    logger.info(f"\t\t{file}\n")

    text = f"Group file {file} correctly read"
    display_alert(text, "success")
    return id_dict, df_group


def index_reference_genome(file, logger):
    suffixes = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    logger.info(f"CHECKING REFERENCE FILE :")
    logger.info(f"\t\tFasta : {file}")
    text = f"Searching index files for reference file {file}"
    display_alert(text, "info")
    cpt_test_files = 0
    for suffix in suffixes:
        tested_file = file + suffix
        if not path.exists(tested_file):
            text = f"Indexation in progress with bwa index..."
            display_alert(text, "secondary")
            cmd = f'bwa index {file}'
            logger.info(f"\t\tbwa index cmd : {cmd}")
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if process.returncode:
                log = f'FAILED EXECUTION : {cmd}\n{process.stdout}\n{process.stderr}'
                raise ValueError(log)
            else:
                # display(msg_button(f"SUCCESSFUL INDEXATION: {cmd}", 'green', 'classic'))
                logger.info(f"\t\tLog : {process.stdout + process.stderr}\n")
            break
        else:
            cpt_test_files += 1  # Avoid to get 1 notification per suffix tested

    if cpt_test_files == len(suffixes):
        text = f"Index files already available ({file}).... Skip Indexation Step"
        at = 'warning'
    else:
        at = 'success'
        text = f"""
### indexation with bwa index Successful 

<hr>

* REFERENCE GENOME : {file}
* COMMAND :  {cmd}
* Log bwa into log file
"""
    display_alert(text, at)


def fastq_stats(fastq_file, stat_dir, logger):
    file = os.path.basename(fastq_file)
    stat_file = os.path.join(stat_dir, file + ".fastqstat")

    if not os.path.exists(stat_file):

        text = f"In progress for {file}..."
        display_alert(text, "secondary")

        cmd = f'fastq-stats -D {fastq_file} > {stat_file}'
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        logger.info(f"\t\t\tfastq-stats cmd : {cmd}")

        if process.returncode:
            text = f"Failed execution ({file}).... see log file, resolve the problem and try again"
            at = 'danger'
        else:
            at = 'success'
            text = f"fastq_stats executed successfully ({stat_file})"
        logger.info(f"\t\t\tLog fastq-stats : {process.stdout + process.stderr}")
    else:
        text = f"Fastq stats previously runned for ({stat_file})"
        at = 'warning'

    display_alert(text, at)


def fastq_stats_dir(fastq_dir, stat_dir, df, logger):
    make_dir(stat_dir)

    logger.info(f"FASTQ STAT :")
    logger.info(f"\t\tfastq stat directory : {stat_dir}")
    text = f"Generating some stat about raw data into ({stat_dir})..."
    display_alert(text, 'info')

    text = f"Reading fastq directory ({fastq_dir}) and group file..."
    display_alert(text, 'info')

    for file_name in os.listdir(fastq_dir):

        if '.fastq' in file_name or '.fq' in file_name or '.fastq.gz' in file_name or '.fq.gz' in file_name:
            read_group = file_name.split("_")[0]

            if read_group in df['sample'].values:
                fastq_stats(os.path.join(fastq_dir, file_name), stat_dir, logger)
            else:
                text = f"Read group {read_group} in this file ({file_name}) not in group file.... Skip fastq stats"
                at = 'warning'
                display_alert(text, at)

    text = f"""
    ### All stat files generated successfully with fastq-stat

    <hr>

    * FASTQ DIR : {fastq_dir}
    * STAT DIR : {stat_dir}
    """

    display_alert(text, "success")


def merge_fastqstat(stat_file, stat_dir, logger):
    if os.path.exists(stat_file):
        os.remove(stat_file)

    # open the file merging the fastq-stat results on the whole set of samples
    f = open(stat_file, "a+")
    header = "sample\tfile\treads\tlen\tlen_mean\tlen_stdev\tlen_min\tphred\twindow-size\tcycle-max\tqual_min\tqual_max\tqual_mean\tqual_stdev\t%A\t%C\t%G\t%T\t%N\ttotal_bases"
    f.write(header + "\n")

    # Parsing the individual fastq-stats file
    for file_name in os.listdir(stat_dir):

        if 'fastqstat' not in file_name:
            continue

        file_base = os.path.basename(file_name).split('.')[0]
        sample = file_base.split('_')[0]
        newLine = f"{sample}\t{file_base}\t"

        with open(os.path.join(stat_dir, file_name), "r") as stat:
            for line in stat:
                line = line.rstrip()
                newLine += line.split('\t')[1] + "\t"
        f.write(newLine + "\n")


def fastq2bam_dir(reference, fastq_dir, df, cpu, bam_dir, logger):
    make_dir(bam_dir)

    logger.info(f"\nMAPPING STEP :")
    logger.info(f"\t\tFastq directory : {fastq_dir}")
    logger.info(f"\t\tMapping directory : {bam_dir}")

    text = f"Starting mapping with bwa mem in {bam_dir}..."
    display_alert(text, 'info')

    for file_name in os.listdir(fastq_dir):

        if '1.fastq' in file_name or '1.fq' in file_name:
            read_group = file_name.split("_")[0]

            if read_group in df['sample'].values:
                fastq2bam(reference, os.path.join(fastq_dir, file_name), cpu, bam_dir, logger)

    text = f"""
    ### Mapping step done with bwa mem and samtools sort

    <hr>

    * FASTQ DIR : {fastq_dir}
    * BAM DIR : {bam_dir}
    """

    display_alert(text, "success")


def fastq2bam(reference, fastq_file, cpu, bam_dir, logger):
    file_name = os.path.basename(fastq_file)
    read_group = file_name.split("_")[0]

    text = f"Mapping in progress for {read_group}..."
    display_alert(text, "secondary")

    fastq2_file = fastq_file.replace("1.f", "2.f")
    sam_file = os.path.join(bam_dir, read_group + ".sam")
    bam_file = os.path.join(bam_dir, read_group + ".bam")

    if not os.path.exists(bam_file):

        cmd1 = f'bwa mem -M -t {cpu} {reference} {fastq_file} {fastq2_file} -o {sam_file}'
        process1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
        logger.info(f"\t\t\tbwa mem cmd : {cmd1}")

        if process1.returncode:
            text = f"Failed execution ({read_group}).... see log file, resolve the problem and try again"
            at = 'danger'
            display_alert(text, at)
            logger.info(f"\t\t\tLog bwa mem : {process1.stdout + process1.stderr}\n")
        else:
            at = 'success'
            text = f"bwa mem executed successfully ({sam_file})"
            display_alert(text, at)
            logger.info(f"\t\t\tLog bwa mem : {process1.stdout + process1.stderr}\n")

        text = f"Sort bam file in progress for {read_group}..."
        display_alert(text, "secondary")

        cmd2 = f'samtools sort {sam_file} -@ {cpu} -o {bam_file} '
        process2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
        logger.info(f"\t\t\tsamtools cmd : {cmd2}")

        if process2.returncode:
            text = f"Failed execution ({sam_file}).... see log file, resolve the problem and try again"
            at = 'danger'
            display_alert(text, at)
        else:
            at = 'success'
            text = f"Sort bam file  successfully executed ({bam_file})"
            display_alert(text, at)
            os.remove(sam_file)
        logger.info(f"\t\t\tLog samtools sort : {process2.stdout + process2.stderr}")

    else:
        text = f"Mapping previously runned for ({bam_file})"
        at = 'warning'
        display_alert(text, at)


def samtools_flagstat(bam_name, stat_dir, logger):
    bam = os.path.basename(bam_name)
    stat_file = os.path.join(stat_dir, bam + ".samtoolsFlagstat")

    text = f"Generating mapping stat (with samtools flagstat) for {bam}..."
    display_alert(text, "secondary")

    if not os.path.exists(stat_file):

        cmd = f'samtools flagstat {bam_name} > {stat_file}'
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        logger.info(f"\t\t\tsamtools cmd : {cmd}")

        if process.returncode:
            text = f"Failed execution ({bam}).... see log file, resolve the problem and try again"
            at = 'danger'
        else:
            at = 'success'
            text = f"samtools flagstat executed successfully ({os.path.basename(stat_file)})"

        logger.info(f"\t\t\tLog samtools : {process.stdout + process.stderr}")
        display_alert(text, at)
    else:
        text = f"Stat previously executed for ({bam})"
        at = 'warning'
        display_alert(text, at)


def bam_to_F0x2_bam_dir(bam_dir, id_dict, cpu, output_dir, logger):
    for id in id_dict:
        bam = os.path.join(bam_dir, id + ".bam")
        new_bam = os.path.join(output_dir, id + "_F0x2.bam")

        if check_file(new_bam) == True:
            text = f"File {new_bam} already existed"
            display_alert(text, 'warning')
        else:
            display_alert(f"Filtering {bam} in progress...", 'secondary')
            cmd = f'samtools view -b -h -F 0x2 -o {new_bam} -@ {cpu} {bam}'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            logger.info(f"\t\t\tsamtools view cmd : {cmd}")

            if process.returncode:
                log = f'FAILED EXECUTION : {cmd}\n{process.stdout}\n{process.stderr}'
                display_alert(log, 'danger')
            else:
                display_alert(f"Filtering bam executed sucessfully", 'success')

            if len(process.stdout) > 0:
                logger.info(f"\t\t\tLog samtools view (STDOUT): {process.stdout}")
            if len(process.stderr) > 0:
                logger.info(f"\t\t\tLog samtools view (STDERR): {process.stderr}")
            if len(process.stdout) == 0 and len(process.stderr) == 0:
                logger.info(f"\t\t\tLog samtools view : ok")

    text = f"""### Extracting unmapped reads
<hr>
* BAM DIR : {bam_dir}
* FILTERED BAM DIR : {output_dir}
"""

    display_alert(text, "success")


def merge_flagstat(output_dir, logger):
    stat_file = os.path.join(output_dir, "all_flagstat.csv")

    if os.path.exists(stat_file):
        os.remove(stat_file)

    logger.info(f"Merging mapping stat into : {stat_file}")
    text = f"Compiling mapping stat into an unique file..."
    display_alert(text, "secondary")

    f = open(stat_file, "w")
    header = "sample,MAPPED,PAIRED,UNMAPPED"
    f.write(header)

    # Parsing the individual fastq-stats file
    for file_name in os.listdir(output_dir):

        if not "samtoolsFlagstat" in file_name:
            continue

        sample = os.path.basename(file_name).split('.')[0]
        newLine = f"\n{sample}"

        with open(os.path.join(output_dir, file_name), "r") as stat:
            for line in stat:
                line = line.rstrip()

                if 'mapped (' in line or 'paired (' in line or 'singleton' in line:
                    # print(line.split('(')[1].split('%')[0])
                    newLine += f",{line.split('(')[1].split('%')[0]}"

        f.write(newLine)

    f.close()
    at = 'success'
    text = f"Compiling successfully done. See {stat_file}"
    display_alert(text, at)


def format_60(txt):
    output = ''
    cpt = 0
    for c in txt:
        if cpt < 60:
            output += c
            cpt += 1
        else:
            output += '\n' + c
            cpt = 1
    return (output)


def check_file(file):
    import os.path
    from os import path
    return (path.isfile(file))


def abyss_pe(project_name, id, k, bam_dir, output_dir, logger):
    tag = project_name + '_' + id + '_' + str(k)
    output_abyss_dir = os.path.join(output_dir, id + "_k" + str(k))
    make_dir(output_abyss_dir)
    test_file = os.path.join(output_abyss_dir, tag + "-contigs.fa")

    if check_file(test_file) == True:
        display_alert(f"File {test_file} already existed", 'warning')
    else:
        bam = os.path.join(bam_dir, id + "_F0x2.bam")
        display_alert(f"Assembly for {bam} ({k}) in progress...", 'secondary')
        cmd = f'abyss-pe -C {output_abyss_dir} name={tag} k={k} in={bam}'
        logger.info(f"\t\t\tABySS cmd : {cmd}")
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if process.returncode:
            log = f'FAILED EXECUTION : {cmd}\n{process.stdout}\n{process.stderr}'
            display_alert(log, 'warning')
        else:
            display_alert(f"Abyss successfully executed", 'success')

        if len(process.stdout) > 0:
            logger.info(f"\t\t\tLog ABySS (STDOUT): {process.stdout}")
        if len(process.stderr) > 0:
            logger.info(f"\t\t\tLog ABySS (STDERR): {process.stderr}")
        if len(process.stdout) == 0 and len(process.stderr) == 0:
            logger.info(f"\t\t\tLog ABySS : ok")


def filter_fastq_threshold(file, output_file, threshold):
    from Bio import SeqIO
    if check_file(output_file) == True:
        display_alert(f"File {output_file} already existed", 'warning')
    else:
        with open(output_file, 'w') as o:
            for seq in SeqIO.parse(file, "fasta"):
                if len(seq) >= threshold:
                    o.write('>' + str(seq.description) + '\n')
                    txt = str(seq.seq)
                    lines = txt.split('\n')
                    for line in lines:
                        if len(line) > 0:
                            if line[0] == '>':
                                o.write(line + '\n')
                            else:
                                o.write(format_60(line) + '\n')


def def_stats():
    stat_len = ["total_length"]
    stats_N_hide = ["mean_length", "longest", "n10", "n20", "n30", "n40"]
    stats_N = ["shortest", "n50", "n60", "n70", "n80", "n90"]
    stats_L_hide = ["n20n", "n30n", "n40n", "n60n", "n70n", "n80n"]
    stats_L = ["n10n", "n50n", "n90n", "number"]
    stats_gap = ["N_count", "Gaps"]
    stats = stat_len + stats_N_hide + stats_N + stats_L_hide + stats_L + stats_gap
    return (stat_len, stats_N_hide, stats_N, stats_L_hide, stats_L, stats_gap, stats)


def create_stats_files(stats, output_dir):
    header = ["id", "k", "stat", "value"]
    for stat in stats:
        output_file = os.path.join(output_dir, "assembly-stats-" + stat + ".csv")
        with open(output_file, 'w') as o:
            o.write('\t'.join(header) + '\n')


def fill_stats_files(input_dir, id, k, output_dir, threshold, logger):
    fasta = id + "_k" + str(k) + "_thr" + str(threshold) + "-contigs.fasta"
    stat_dict = parse_assembly_stats_adapted(os.path.join(input_dir, fasta), logger)
    for stat in stat_dict:
        with open(os.path.join(output_dir, "assembly-stats-" + stat + ".csv"), 'a') as o:
            o.write('\t'.join([id, str(k)]) + '\t' + '\t'.join([stat, stat_dict[stat]]) + '\n')


def parse_assembly_stats_adapted(file, logger):
    cmd = f'assembly-stats -s {file}'
    text = f"Generating stat (with assembly-stats) for {file}..."
    display_alert(text, "secondary")
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    logger.info(f"\t\t\tAssembly-stats cmd : {cmd}")

    if process.returncode:
        log = f'FAILED EXECUTION : {cmd}\n{process.stdout}\n{process.stderr}'
        display_alert(log, 'warning')
    else:
        display_alert(f"assembly-stats successfully done ({file})", 'success')
    if len(process.stdout) > 0:
        logger.info(f"\t\t\tLog assembly-stats (STDOUT): {process.stdout}")
    if len(process.stderr) > 0:
        logger.info(f"\t\t\tLog assembly-stats (STDERR): {process.stderr}")
    if len(process.stdout) == 0 and len(process.stderr) == 0:
        logger.info(f"\t\t\tLog assembly-stats : ok")

    input = process.stdout.split('\n')
    stats = {}
    for line in input:
        if len(line) > 0:
            fields = line.split('\t')  # [0] fasta, [1] stat, [2] val
            stats[fields[1]] = fields[2]
    return stats

def fasta_stats(fasta_dir, stat_file, logger):

    text = f"Generating abyss stat (with assembly-stat) for fasta files in {fasta_dir}..."
    display_alert(text, "secondary")

    if not os.path.exists(stat_file):

        cmd = f'assembly-stats -t {fasta_dir}/*fasta > {stat_file}'
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        logger.info(f"\t\t\tassembly-stats cmd : {cmd}")

        if process.returncode:
            text = f"Failed execution.... see log file, resolve the problem and try again"
            at = 'danger'
        else:
            at = 'success'
            text = f"stat executed successfully"

        logger.info(f"\t\t\tLog assembly-stats : {process.stdout + process.stderr}")
        display_alert(text, at)
    else:
        text = f"Stat previously executed with ({stat_file})"
        at = 'warning'
        display_alert(text, at)

# def parse_assembly_stats(file) :
#     input = !assembly-stats -s $file
#     stats = []
#     for line in input : 
#         fields = line.split('\t') # [0] fasta, [1] stat, [2] val
#         stats.append(fields[2])
#     return(stats)

def index_blast(db_file,logger):

    suffixes = [".nhr", ".nin", ".nsq"]
    logger.info(f"CHECKING BLAST INDEX FOR VECTOR BANK :")
    logger.info(f"\t\tFasta : {db_file}")
    text = f"Searching index files for reference file {db_file}"
    display_alert(text, "info")
    cpt_test_files = 0

    for suffix in suffixes:
        tested_file = db_file + suffix

        if not path.exists(tested_file):
            text = f"Indexation in progress with makeblastdb..."
            display_alert(text, "secondary")
            cmd = f"makeblastdb -in {db_file} -dbtype 'nucl' -blastdb_version 4"
            logger.info(f"\t\tmakeblastdb cmd : {cmd}")
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if process.returncode:
                log = f'FAILED EXECUTION : {cmd}\n{process.stdout}\n{process.stderr}'
                raise ValueError(log)
            else:
                # display(msg_button(f"SUCCESSFUL INDEXATION: {cmd}", 'green', 'classic'))
                logger.info(f"\t\tLog : {process.stdout + process.stderr}\n")
            break
        else:
            cpt_test_files += 1  # Avoid to get 1 notification per suffix tested

    if cpt_test_files == len(suffixes):
        text = f"Index files already available ({db_file}).... Skip Indexation Step"
        at = 'warning'
    else:
        at = 'success'
        text = f"""
### indexation with blast db index Successful 

<hr>

* BANK : {db_file}
* COMMAND :  {cmd}
* Log blast into log file
"""
    display_alert(text, at)


def run_vecscreen(fasta_dir, id, k, threshold, output_dir, bank,logger):
    from Bio import SeqIO

    # print("--- id : " + id)
    tag = id + "_k" + str(k) + "_thr" + str(threshold)
    fasta_file = os.path.join(fasta_dir, tag + "-contigs.fasta")
    output_file = os.path.join(output_dir, tag + '_vecscreen.fasta')

    if path.exists(fasta_file) and not path.exists(output_file):

        nb_seq, nb_seq_no_hits, nb_seq_suspicious = 0, 0, 0
        text = f"VecScreen for {id} (k = {k}) in progress ..."
        display_alert(text, "secondary")

        with open(output_file, "w") as o:

            for seq in SeqIO.parse(fasta_file, "fasta"):
                # Vecscreen only works with fasta file, not sequence
                tmp_file = output_dir + "tmp.fasta"

                with open(tmp_file, "w") as tmp:
                    tmp.write('>' + str(seq.description) + '\n')
                    tmp.write(format_60(str(seq.seq)) + '\n')

                cmd = f'vecscreen -d {bank} -i {tmp_file} -f 3'
                process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                logger.info(f"\t\t\tvecscreen cmd : {cmd}")

                if process.returncode:
                    text = f"Failed execution ({cmd}).... see log file, resolve the problem and try again"
                    at = 'danger'
                else:
                    at = 'success'
                    text = f"vecscreen executed successfully ({os.path.basename(tmp_file)})"

                    if 'No hits found' in process.stdout:
                        o.write('>' + str(seq.description) + '\n')
                        o.write(format_60(str(seq.seq)) + '\n')
                        nb_seq += 1
                        nb_seq_no_hits += 1
                    else:
                        # print(process.stdout)
                        # print("Hits : >" + str(seq.description))
                        nb_seq_suspicious += 1
                        bounds = []
                        for line in process.stdout.split("\n"):
                            if len(line) > 0:
                                if line[0].isnumeric():
                                    for value in line.split('\t'):
                                        bounds.append(int(value))
                        min_seq = 0
                        max_seq = len(seq)
                        min_bounds = min(bounds)
                        max_bounds = max(bounds)
                        if min_bounds == min_seq + 1:
                            min_seq = max_bounds
                            new_length = max_seq - max_bounds
                            if new_length >= threshold:
                                header_fields = str(seq.description).split()
                                o.write('>' + header_fields[0] + " " + str(new_length) + " " + " ".join(
                                    header_fields[2:]) + '\n')
                                o.write(format_60(str(seq.seq[min_seq:max_seq])) + '\n')
                                nb_seq += 1
                        elif max_bounds == max_seq:
                            max_seq = min_bounds
                            new_length = max_seq
                            if new_length >= threshold:
                                header_fields = str(seq.description).split()
                                o.write('>' + header_fields[0] + " " + str(new_length) + " " + " ".join(
                                    header_fields[2:]) + '\n')
                                o.write(format_60(str(seq.seq[min_seq:max_seq])) + '\n')
                                nb_seq += 1
                        else:
                            print("Warning (seq : " + str(
                                seq.description) + ") : the detected contamination appears to be in the middle of the contig !")
                    cmd = f'rm {tmp_file}'
                    # process = f'subprocess.run(cmd, shell=True, capture_output=True, text=True)'

        text = f"""
        ### VecScreen
        <hr>
        * NO HITS : {nb_seq_no_hits}index_blast(vec_file.value)
        * SUSPICIOUS : {nb_seq_suspicious}
        * KEPT : {nb_seq}
        """
        display_alert(text, "success")


def copy_cluster(dir_ctg, dir_clustering):
    text = f"Copying all contigs from each sample into dir_clustering on progress..."
    display_alert(text, "secondary")

    for filename in os.listdir(dir_ctg):

        # output_sample = os.path.join(dir_ctg, filename)
        readgroup = filename.split('_')[0]

        target = os.path.join(dir_clustering, filename)
        shutil.copyfile(os.path.join(dir_ctg, filename), target)

        sedcmd = f"sed 's/>/>{readgroup}_/' {target} -i"
        process = subprocess.run(sedcmd, shell=True, capture_output=True, text=True)
        # logger.info(f"\t\t\tsed cmd : {sedcmd}")

        if process.returncode:
            text = f"Renaming Failed.... see log file, resolve the problem and try again"
            at = 'danger'
            display_alert(text, at)

    text = f"Copying all contigs into dir_clustering done"
    display_alert(text, "info")


def merging_cluster(dir_clustering, contigs_file):
    if os.path.exists(contigs_file):
        os.remove(contigs_file)

    merged = open(contigs_file, 'w')

    text = f"Merging all contigs from each sample into one file... on progress"
    display_alert(text, "secondary")

    for filename in os.listdir(dir_clustering):

        if "_allContigs.fa" in filename:
            # print("...."+filename+"\n")
            continue

        # print("----"+filename+"\n")
        fasta = open(os.path.join(dir_clustering, filename), 'r')
        merged.write(fasta.read())
        fasta.close()
        os.remove(os.path.join(dir_clustering, filename))

    merged.close()

    text = f"fasta file successfully done {contigs_file}"
    display_alert(text, "success")

def cdhit(contigs, c, s, cdhit_file, logger):

    text = f"Clustering all contigs with cdhit-est.. on progress"
    display_alert(text, "secondary")

    if check_file(cdhit_file) == True:
        display_alert(f"File {cdhit_file} already existed", 'warning')
    else:
        cmd = f"cd-hit-est -c {c} -s {s} -i {contigs} -o {cdhit_file}  "

        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        logger.info(f"\t\t\tcdhit cmd : {cmd}")

        if process.returncode:
            text = f"Clustering Failed.... see log file, resolve the problem and try again"
            # print({process.stdout + process.stderr})
            at = 'danger'
        else:
            text = f"Clustering successfully done !"
            at = 'success'
        display_alert(text, at)
        logger.info(f"\t\t\tcdhit log : {process.stdout + process.stderr}")


def parse_cdhit(cdhit_cluster, df_group, cdhit_csv):
    cluster = {}

    sp = None
    tag = None
    ctg = None
    sp = None
    sp_list = None
    ctg_list = None

    text = f"Merging clustering information for {cdhit_cluster}..."
    display_alert(text, "secondary")

    with open(cdhit_cluster, "r") as stat:

        for line in stat:

            if ">Cluster" in line:

                # With the precedent cluster
                if tag is not None:
                    cluster[tag] = {
                        'pb': ln,
                        'sp': sp,
                        'ctg_list': ctg_list[:-1],
                        'sp_list': sp_list[:-1]
                    }

                # Initialization of variables for the new cluster
                sp = None
                tag = None
                ctg = None
                sp = None
                sp_list = ""
                ctg_list = ""
                ln = None

            else:
                ctg = line.split(">")[1].split(".")[0]
                sample = ctg.split("_")[0]

                ctg_list += ctg + ","
                sp_tmp = ""

                if (df_group[(df_group['sample'] == sample)]['Species'].count() == 1):
                    sp_tmp = [value for key, value in df_group[(df_group['sample'] == 'AB')]['Species'].items()][0]

                sp_list += sp_tmp + ","

                if "*" in line and "nt" in line:
                    ln = line.split("n")[0].split("\t")[1]
                    tag = ctg
                    sp = sp_tmp
                    # print(f"TAG : {tag} {sp}")

    cluster[tag] = {
        'pb': ln,
        'sp': sp,
        'ctg_list': ctg_list[:-1],
        'sp_list': sp_list[:-1]
    }

    # Save into  csv file
    # Create dataframe from dic and make keys, index in dataframe
    df_cluster = pd.DataFrame.from_dict(cluster, orient='index', columns=['pb', 'sp', 'ctg_list', 'sp_list'])
    df_cluster.reset_index(level=0, inplace=True)
    df_cluster['ln'] = (df_cluster['sp_list'].str.count(',') + 1)
    df_cluster.to_csv(cdhit_csv, index=False)

    at = 'success'
    text = f"Compiling successfully done. See {cdhit_csv}"
    display_alert(text, at)

    return df_cluster


def samtools_index(bam_name, logger):
    bam = os.path.basename(bam_name)

    text = f"Bam indexing {bam}..."
    display_alert(text, "secondary")

    cmd = f'samtools index {bam_name}'
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    logger.info(f"\t\t\tsamtools cmd : {cmd}")

    if process.returncode:
        text = f"Failed execution ({bam}).... see log file, resolve the problem and try again"
        at = 'danger'
    else:
        at = 'success'
        text = f"samtools index executed successfully"
    logger.info(f"\t\t\tLog samtools : {process.stdout + process.stderr}")
    display_alert(text, at)


def get_seq_nb(fasta):

    result = subprocess.run(['grep', '-c', '>', fasta], capture_output=True, text=True)
    return (result.stdout)


def generating_panref(ref_fasta, ctg_fasta, panref_fasta, logger):

    filenames = [ref_fasta, ctg_fasta]

    with open(panref_fasta, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    logger.info(f"\t\t\tgenerating panref fasta : {panref_fasta}")
    display_alert(f"Panreference file successfully created : {panref_fasta}", "success")
    index_reference_genome(panref_fasta, logger)

def anchoring(output_panrefmapping_dir, panrefposi_file, depth, logger):

    logger.info(f"ANCHORING OF CONTIGS ON REFERENCE :")
    logger.info(f"\t\tWorkig directory : {output_panrefmapping_dir}")

    text = f"Anchoring indexing..."
    display_alert(text, "secondary")

    cmd = f'tools/parseBamv7.py -b {output_panrefmapping_dir} -d {depth} -o {panrefposi_file}'
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    logger.info(f"\t\t\placement cmd : {cmd}")

    if process.returncode:
        text = f"Failed execution.... see log file, resolve the problem and try again"
        at = 'danger'
    else:
        at = 'success'
        text = f"script executed successfully"
    logger.info(f"\t\t\tLog : {process.stdout + process.stderr}")
    display_alert(text, at)


def parse_anchoring(ctg_fasta, panrefposi_file, pass_file, nopass_file, length, logger):

    df_anc = pd.read_csv(panrefposi_file, index_col=False, sep=";")
    stat_dict = { }

    # Columns To print
    col_toprint = ['CTG_name', 'CTG_len', 'CHR_name', 'ANCHORING_TAG', 'SAMPLE', 'SAMPLE_list', 'TOT_Hit', 'START_Hit',
                   'START_CHR_posi_sd', 'START_CHR_Min', 'END_Hit', 'END_CHR_posi_sd', 'END_CHR_Max', 'Chr_Ln']

    # Contigs number total
    #nb_tot_ctgs = int(get_seq_nb(ctg_fasta))

    ######################### STATS ON CTGS PLACED JUST ON A SINGLE CHR
    ### ALL CTGS PLACED (with redondant ctgs)
    ctgs_anchored = df_anc.CTG_name.nunique(dropna=True)
    stat_dict['contigs_anchored'] = ctgs_anchored

    #print("------------ Number of Contigs anchored : ", ctgs_anchored, ",", round(ctgs_anchored / nb_tot_ctgs * 100, 2),"% of all ctgs")

    ### JUST ATGS PLACE ON A SINGLE CHR
    # remove toutes les rows avec plusieurs fois le meme contigs
    df_uniq_anc = df_anc.drop_duplicates('CTG_name', keep=False)

    # STATS ON CTGS PLACED JUST ONCE ON A SINGLE CHROMOSOME
    oneside_ctgs_nb = len(df_uniq_anc[df_uniq_anc.ANCHORING_TAG == "5\'"]) + len(
        df_uniq_anc[df_uniq_anc.ANCHORING_TAG == "3\'"])
    botside_Ctgs_nb = len(df_uniq_anc[df_uniq_anc.ANCHORING_TAG == "5\'&3\'"])

    total_ctgs_nb = oneside_ctgs_nb + botside_Ctgs_nb
    stat_dict['uniq_anchored'] = total_ctgs_nb
    stat_dict['uniq_anchored_mean'] = int(df_uniq_anc["CTG_len"].mean())
    stat_dict['uniq_anchored_median'] = int(df_uniq_anc["CTG_len"].median())

    # print("\n------------ Number of contigs anchored at mulitple position :", ctgs_anchored - total_ctgs_nb)
    # print("\n------------ Number of contigs anchored at an unique position : ", total_ctgs_nb, ",",
    #      round(total_ctgs_nb / nb_tot_ctgs * 100, 2), "% of all ctgs")
    # print(f"mean  & median ctg length :", int(df_uniq_anc["CTG_len"].mean()), int(df_uniq_anc["CTG_len"].median()))

    ######################### ADD CHR SD AS COLUMN
    df_uniq_anc['START_CHR_sd'] = df_uniq_anc.START_CHR_posi_sd.str.extract('(.+)\(.*').fillna(0).astype(float)
    df_uniq_anc['END_CHR_sd'] = df_uniq_anc.END_CHR_posi_sd.str.extract('(.+)\(.*').fillna(0).astype(float)

    #### ################## FILTER CTGS PLACED UNIQUELY ON ONE CHR
    filtered_ctgs_nb = 0

    #################### FILTERING Ctgs 5&3
    both_bool = (df_uniq_anc.ANCHORING_TAG == "5\'&3\'")
    nan_bool1 = (df_uniq_anc.START_CHR_sd == -9999.0)
    nan_bool2 = (df_uniq_anc.END_CHR_sd == -9999.0)
    filtered_ctgs_nb = df_uniq_anc[(both_bool) & ((df_uniq_anc.START_CHR_sd < length & ~nan_bool1) | (
                (df_uniq_anc.END_CHR_sd < length) & ~nan_bool2))].CTG_name.count()

    ################### FILTERINg 5 or 3
    one5_bool = (df_uniq_anc.ANCHORING_TAG == "5\'")
    one3_bool = (df_uniq_anc.ANCHORING_TAG == "3\'")

    filtered_ctgs_nb += df_uniq_anc[one3_bool & ~(nan_bool2) & (df_uniq_anc.END_CHR_sd < length)].CTG_name.count() + \
                        df_uniq_anc[one5_bool & ~(nan_bool1) & (df_uniq_anc.START_CHR_sd < length)].CTG_name.count()

    stat_dict['uniq_pass'] = filtered_ctgs_nb

    #print("\n------------ FILTERED CONTIGS : ", filtered_ctgs_nb, ", ", round(filtered_ctgs_nb / nb_tot_ctgs * 100, 2),
          "% of all contigs")

    #################### SAVE FILTERED CONTIGS INTO A FILE
    df_thebest = pd.concat([df_uniq_anc[one3_bool & ~(nan_bool2) & (df_uniq_anc.END_CHR_sd < length)][col_toprint],
                            df_uniq_anc[one5_bool & ~(nan_bool1) & (df_uniq_anc.START_CHR_sd < length)][col_toprint],
                            df_uniq_anc[(both_bool) & ((df_uniq_anc.START_CHR_sd < length & ~nan_bool1) | (
                                        (df_uniq_anc.END_CHR_sd < length) & ~nan_bool2))][col_toprint]])
    df_thebest.to_csv(pass_file, index=False)
    #print(f"Min - Max : {df_thebest.CTG_len.min()} - {df_thebest.CTG_len.max()} \nMean, median length : {int(df_thebest.CTG_len.mean())} - {int(df_thebest.CTG_len.median())}")

    ctg_list = df_thebest.CTG_name
    df_thebad = df_anc[~(df_anc.CTG_name.isin(ctg_list))][col_toprint]
    df_thebad.to_csv(nopass_file, index=False)
    #print(pass_file,nopass_file)

    return stat_dict