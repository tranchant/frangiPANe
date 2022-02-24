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
        text = f"""### {at }
<hr>
Directory {path} already existed"""

    display_alert(text,at)


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
    text=f"Searching index files for reference file {file}"
    display_alert(text, "info")
    cpt_test_files = 0
    for suffix in suffixes:
        tested_file = file + suffix
        if not path.exists(tested_file):
            text=f"Indexation in progress with bwa index..."
            display_alert(text, "secondary")
            cmd = f'bwa index {file}'
            logger.info(f"\t\tbwa index cmd : {cmd}")
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if process.returncode:
                log = f'FAILED EXECUTION : {cmd}\n{process.stdout}\n{process.stderr}'
                raise ValueError(log)
            else:
                #display(msg_button(f"SUCCESSFUL INDEXATION: {cmd}", 'green', 'classic'))
                logger.info(f"\t\tLog : {process.stdout + process.stderr}\n")
            break
        else :
            cpt_test_files += 1 # Avoid to get 1 notification per suffix tested

    if cpt_test_files == len(suffixes) :
        text = f"Index files already available ({file}).... Skip Indexation Step"
        at='warning'
    else:
        at='success'
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
    stat_file = stat_dir + file + ".fastqstat"

    if not os.path.exists(stat_file):

        text = f"In progress for {file}..."
        display_alert(text, "secondary")

        cmd = f'fastq-stats -D {fastq_file} > {stat_file}'
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        logger.info(f"\t\t\tfastq-stats cmd : {cmd}")

        if process.returncode:
            text = f"Failed execution ({file}).... see log file, resolve the problem and try again"
            at='danger'
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


def merge_fastqstat(stat_file, stat_dir,logger):

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

        with open(os.path.join(stat_dir,file_name), "r") as stat:
            for line in stat:
                line = line.rstrip()
                newLine += line.split('\t')[1] + "\t"
                #if 'total' in line.split('\t')[0]:
                    #if sample not in total_base:
                    #    total_base[sample] = 0
                    #total_base[sample] += int(line.split('\t')[1])
        f.write(newLine + "\n")

    #total_base_value = list(total_base.values())
    #sample = [*total_base]
    #return total_base_value, sample


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
                fastq2bam(reference, os.path.join(fastq_dir, file_name), cpu,bam_dir, logger)

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

    text=f"Mapping in progress for {read_group}..."
    display_alert(text, "secondary")

    fastq2_file = fastq_file.replace("1.f", "2.f")
    sam_file = bam_dir + read_group + ".sam"
    bam_file = bam_dir + read_group + ".bam"

    if not os.path.exists(bam_file) :

        cmd1 = f'bwa mem -M -t {cpu} {reference} {fastq_file} {fastq2_file} -o {sam_file}'
        process1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
        logger.info(f"\t\t\tbwa mem cmd : {cmd1}")

        if process1.returncode:
            text = f"Failed execution ({read_group}).... see log file, resolve the problem and try again"
            at='danger'
            display_alert(text, at)
            logger.info(f"\t\t\tLog bwa mem : {process1.stdout + process1.stderr}\n")
        else:
            at = 'success'
            text = f"bwa mem executed successfully ({sam_file})"
            display_alert(text, at)
            logger.info(f"\t\t\tLog bwa mem : {process1.stdout + process1.stderr}\n")

        #display(msg_button(f"MAPPING STEP FOR ({read_group})... samtools sort in progress", 'blue', 'classic'))
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


#def fastq_to_bam(reference_genome, fastq_dir, id, cpu, output_dir,logger):

#    logger.info(f"\t\tread group: {id}")
#    sam_file = output_dir+id+".sam"
#    bam_file = output_dir+id+".bam"

#    display(f"MAPPING STEP FOR ({id})... bwa mem in progress",'secondary'))
#    cmd1= f'bwa mem -p -M -t { cpu } { reference_genome}  { fastq_dir+"/"+id+".fastq" } -o { sam_file}'
#    process1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
#    logger.info(f"\t\t\tbwa mem cmd : {cmd1}")

#    if process1.returncode:
#        text = f'FAILED EXECUTION : {cmd1} {process1.stdout} { process1.stderr}'
#        at = 'danger'
#        display(text,at)
#        logger.info(f"\t\t\tLog bwa mem : {process1.stdout + process1.stderr}")
#    else:
#        display(f"SUCCESSFUL MAPPING: {id}\n {cmd1}", 'success'))
#        logger.info(f"\t\t\tLog bwa mem : {process1.stdout + process1.stderr}")

#        display(msg_button(f"MAPPING STEP FOR ({id})... samtools sort in progress", 'blue', 'classic'))
#        cmd2 = f'samtools sort {sam_file} -@ {cpu} -o {bam_file} '
#        process2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
#        logger.info(f"\t\t\tsamtools cmd : {cmd2}")

#       if process2.returncode:
#            log = f'FAILED EXECUTION : {cmd2} {process2.stdout} {process2.stderr}'
#            display(msg_button(log, 'warning', 'warning'))
#        else:
#            display(msg_button(f"SUCCESSFUL SAM TO BAM: {id} {cmd2}", 'green', 'classic'))
#            os.remove(sam_file)
#        logger.info(f"\t\t\tLog bwa mem : {process2.stdout + process2.stderr}")


def samtools_flagstat(bam_name, stat_dir, logger):

    bam = os.path.basename(bam_name)
    stat_file = stat_dir + bam + ".samtoolsFlagstat"

    text=f"Generating mapping stat (with samtools flagstat) for {bam}..."
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
        bam = bam_dir + id + ".bam"
        new_bam = output_dir + id + "_F0x2.bam"

        if check_file(new_bam) == True :
            text=f"File {new_bam} already existed"
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

            if len(process.stdout) > 0 :
                logger.info(f"\t\t\tLog samtools view (STDOUT): {process.stdout}")
            if len(process.stderr) > 0 :
                logger.info(f"\t\t\tLog samtools view (STDERR): {process.stderr}")
            if len(process.stdout) == 0 and len(process.stderr) == 0 :
                logger.info(f"\t\t\tLog samtools view : ok")


    text = f"""### Extracting unmapped reads
<hr>
* BAM DIR : {bam_dir}
* FILTERED BAM DIR : {output_dir}
"""

    display_alert(text, "success")


def merge_flagstat(output_dir, logger):

    stat_file = output_dir + "all_flagstat.csv"

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
    for c in txt :
        if cpt < 60 :
            output += c
            cpt += 1
        else :
            output += '\n' + c
            cpt = 1
    return(output)


def check_file(file):
    import os.path
    from os import path
    return(path.isfile(file))


def abyss_pe(project_name, id, k, bam_dir, output_dir, logger):
    output_abyss_dir = os.path.join(output_dir, id + "_k" + str(k) )
    make_dir(output_abyss_dir)
    test_file = output_abyss_dir + project_name + '_' + id + '_' + str(k) + "-contigs.fa"

    if check_file(test_file) == True:
        display_alert(f"File {test_file} already existed", 'warning')
    else:
        bam = id + "_F0x2.bam"
        display_alert(f"Assembly for {bam} ({k}) in progress...",'secondary')
        cmd = f'abyss-pe -C {output_abyss_dir} name={project_name}_{id}_{str(k)} k={k} in={bam_dir}{bam}'
        logger.info(f"\t\t\tABySS cmd : {cmd}")
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if process.returncode:
            log = f'FAILED EXECUTION : {cmd}\n{process.stdout}\n{process.stderr}'
            display_alert(log, 'warning')
        else:
            display_alert(f"Abyss successfully executed",'success')

        if len(process.stdout) > 0:
            logger.info(f"\t\t\tLog ABySS (STDOUT): {process.stdout}")
        if len(process.stderr) > 0:
            logger.info(f"\t\t\tLog ABySS (STDERR): {process.stderr}")
        if len(process.stdout) == 0 and len(process.stderr) == 0:
            logger.info(f"\t\t\tLog ABySS : ok")


def filter_fastq_threshold(file, output_file, threshold):

    from Bio import SeqIO
    if check_file(output_file) == True :
        display_alert(f"File {output_file} already existed",  'warning')
    else :
        with open(output_file, 'w') as o :
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
    stats_L_hide = ["n20n",  "n30n", "n40n", "n60n", "n70n", "n80n"]
    stats_L = ["n10n", "n50n", "n90n", "number"]
    stats_gap = ["N_count", "Gaps"]
    stats = stat_len + stats_N_hide + stats_N + stats_L_hide + stats_L + stats_gap
    return(stat_len, stats_N_hide, stats_N, stats_L_hide, stats_L, stats_gap, stats)


def create_stats_files(stats, output_dir) :
    header = ["id", "k", "stat", "value"]
    for stat in stats :
        output_file = output_dir + "assembly-stats-" + stat + ".csv"
        with open(output_file, 'w') as o :
            o.write('\t'.join(header) + '\n')


def fill_stats_files(input_dir, id, k, output_dir, threshold, logger) :
    stat_dict = parse_assembly_stats_adapted(input_dir + id + "_k" + str(k) + "_thr" + str(threshold) + ".fasta", logger)
    for stat in stat_dict :
        with open(output_dir + "assembly-stats-" + stat + ".csv", 'a') as o :
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
    if len(process.stdout) > 0 :
        logger.info(f"\t\t\tLog assembly-stats (STDOUT): {process.stdout}")
    if len(process.stderr) > 0 :
        logger.info(f"\t\t\tLog assembly-stats (STDERR): {process.stderr}")
    if len(process.stdout) == 0 and len(process.stderr) == 0 :
        logger.info(f"\t\t\tLog assembly-stats : ok")

    input = process.stdout.split('\n')
    stats = {}
    for line in input :
        if len(line) > 0 :
            fields = line.split('\t') # [0] fasta, [1] stat, [2] val
            stats[fields[1]] = fields[2]
    return stats

# def parse_assembly_stats(file) :
#     input = !assembly-stats -s $file
#     stats = []
#     for line in input : 
#         fields = line.split('\t') # [0] fasta, [1] stat, [2] val
#         stats.append(fields[2])
#     return(stats)


def copy_cluster(dir_abyss_ctg, dir_clustering):

    text = f"Copying all contigs from each sample into dir_clustering on progress..."
    display_alert(text, "secondary")

    for dirname in os.listdir(dir_abyss_ctg):

        if os.path.isdir(dirname):
            continue

        output_sample = os.path.join(dir_abyss_ctg, dirname)
        readgroup = dirname.split('_')[0]
        for filename in os.listdir(output_sample):

            if not "-contigs.fa" in filename:
                continue

            target = os.path.join(dir_clustering, filename)
            shutil.copyfile(os.path.join(output_sample, filename), target)

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

    merged = open(contigs_file, 'a+')

    text = f"Merging all contigs from each sample into one file... on progress"
    display_alert(text, "secondary")

    for filename in os.listdir(dir_clustering):

        if not "-contigs.fa" in filename:
            continue

        fasta = open(os.path.join(dir_clustering, filename), 'r')
        merged.write(fasta.read())
        fasta.close()
        os.remove(os.path.join(dir_clustering, filename))

    merged.close()

    text = f"fasta file succesfully done {contigs_file}"
    display_alert(text, "success")


def cdhit(contigs, c, s, cdhit_file, logger):
    text = f"Clustering all contigs with cdhit-est.. on progress"
    display_alert(text, "secondary")

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
    logger.info(f"\t\t\tscdhit log : {process.stdout + process.stderr}")


def parse_cdhit(cdhit_cluster, df_group, cdhit_csv):
    cluster = {}

    sp = None
    tag = None
    ctg = None
    sp = None
    sp_list = None
    ctg_list = None

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

    # Save into  csv file
    # Create dataframe from dic and make keys, index in dataframe
    df_cluster = pd.DataFrame.from_dict(cluster, orient='index', columns=['pb', 'sp', 'ctg_list', 'sp_list'])
    df_cluster.reset_index(level=0, inplace=True)
    df_cluster['ln'] = (df_cluster['sp_list'].str.count(',') + 1)
    df_cluster.to_csv(cdhit_csv, index=False)

    return df_cluster

def samtools_index(bam_name, logger):

    bam = os.path.basename(bam_name)

    text=f"Bam indexing {bam}..."
    display_alert(text, "secondary")

    cmd = f'samtools index {bam_name}'
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    logger.info(f"\t\t\tsamtools cmd : {cmd}")

    if process.returncode:
        text = f"Failed execution ({bam}).... see log file, resolve the problem and try again"
        at = 'danger'
    else:
        at = 'success'
        text = f"samtools index executed successfully)"
    logger.info(f"\t\t\tLog samtools : {process.stdout + process.stderr}")
    display_alert(text, at)