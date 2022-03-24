#!/usr/bin/env python3

# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html) and GPLv3
# Intellectual property belongs to IRD, ... and SouthGreen development platform
# Written by Clothilde Chenal and Christine Tranchant-Dubreuil
# Copyright 2021

"""

    Use it to import functions to run bioinformatics softwares.

    :Example:

    >>> from tools import jupyter

"""

import ipywidgets as widgets

from IPython.display import display, HTML
from IPython.core.magic import register_cell_magic

import pickle
import pandas as pd
import panel as pn
import logging, os
from time import localtime, strftime

import param
import matplotlib.pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import seaborn as sns

import subprocess

# save variables into a file
file_load = 'frangiPANe.p'


@register_cell_magic
def bgc(color, cell=None):
    script = (
        "var cell = this.closest('.jp-CodeCell');"
        "var editor = cell.querySelector('.jp-Editor');"
        "editor.style.background='{}';"
        "this.parentNode.removeChild(this)"
    ).format(color)

    display(HTML('<img src onerror="{}">'.format(script)))


def init_log(output_dir, project_name):
    start_time = strftime("%d-%m-%Y_%H:%M:%S", localtime())

    log_file = f'frangiPANe_{project_name}-{start_time}.log'
    log_file_path = os.path.join(output_dir, log_file)

    # Create a custom logger
    logger = logging.getLogger()

    # set logger level at DEBUG to print all
    logger.setLevel(logging.DEBUG)

    # Create formatters and add it to handlers
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s', "%Y-%m-%d %H:%M:%S")

    # to log messages to a file
    fh = logging.FileHandler(log_file_path, mode='w')
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    at = "success"
    text = f"Log {log_file_path} created"
    display_alert(text, at)
    return logger


def add_variable(file2save, key, value):
    # get variables stored into the file
    file = open(file2save, 'rb')
    variables = pickle.load(file)
    file.close()
    variables[key] = value

    # add variables and save into the file
    file2 = open(file2save, 'wb')
    pickle.dump(variables, file2)
    file2.close()


def print_variables(file2read):
    file = open(file2read, 'rb')
    variables = pickle.load(file)
    print(variables)
    file.close()


def get_variables(file2read, key):
    file = open(file2read, 'rb')
    variables = pickle.load(file)
    file.close()

    value = 0
    if key in variables:
        value = variables[key]
    return value


def add_css():
    from IPython.core.display import HTML
    css_warning = "<style>.mywarning { color:white ; font-size:100%; font-weight: bold;}"
    HTML("""
    <style>
    """ + css_warning + """
    </style>
    """)
    css_classic = "<style>.myclassic { color:#5B5B5B; font-size:100%; font-weight: bold;}"
    HTML("""
    <style>
    """ + css_classic + """
    </style>
    """)


def widget_display(df):
    out = widgets.Output()
    with out:
        display(df)
    return out


def display_alert(text, at):
    display(pn.pane.Alert(text.format(alert_type=at), alert_type=at))


def box_config():
    # cmd result
    text = "No filled"
    at = 'warning'
    result = pn.pane.Alert(text.format(alert_type=at), alert_type=at, height=200)

    # button
    print_btn = pn.widgets.Button(name='SAVE', width=100, button_type='primary')
    init_btn = pn.widgets.Button(name='INIT', width=100, button_type='primary')
    load_btn = pn.widgets.Button(name='LOAD', width=100, button_type='primary')

    # form
    project_namef = pn.widgets.TextInput(name='Project name :', placeholder='Enter the name here...')
    out_dirf = pn.widgets.TextInput(name='Output directory :', placeholder='Enter the directory path here...')
    fastq_dirf = pn.widgets.TextInput(name='Fastq directory', placeholder='Enter the directory path here...')
    group_filef = pn.widgets.TextInput(name='Group file', placeholder='Enter the file path here...')
    ref_filef = pn.widgets.TextInput(name='Reference file', placeholder='Enter the file path here...')
    vec_filef = pn.widgets.TextInput(name='Univec file', placeholder='Enter the file path here...')
    cpu_pf = pn.widgets.IntSlider(name='CPU number', start=1, end=36, step=1, value=6)

    def reinit_form(event):
        project_namef.value = ""
        out_dirf.value = ""
        fastq_dirf.value = ""
        group_filef.value = ""
        ref_filef.value = ""
        vec_filef.value = ""
        cpu_pf.value = 6

        at = 'danger'
        text = f"""
        ### WARNING : Fill form
        """

        text = "empty"

        return

    def load_form(event):

        file = open(file_load, 'rb')
        variables = pickle.load(file)
        project_namef.value = variables["project_name"]
        out_dirf.value = variables["out_dir"]
        fastq_dirf.value = variables["fastq_dir"]
        group_filef.value = variables["group_file"]
        ref_filef.value = variables["ref_file"]
        vec_filef.value = variables["vec_file"]
        cpu_pf.value = variables["cpu"]
        file.close()

        text = "empty"
        return

    def print_value(event):
        project_name = project_namef.value.rstrip('\s')
        output_dir = out_dirf.value.rstrip('\s')
        fastq_dir = fastq_dirf.value.rstrip('\s')
        group_file = group_filef.value.rstrip('\s')
        ref_file = ref_filef.value.rstrip('\s')
        vec_file = vec_filef.value.rstrip('\s')
        cpu = cpu_pf.value

        if not output_dir or not project_name or not fastq_dir or not group_file or not ref_file or not vec_file:
            at = 'danger'
            text = f"""
    ### WARNING : Fields empty !
    
    * PROJECT NAME : {project_name}
    * OUTPUT DIRECTORY : {output_dir}
    * FASTQ DIRECTORY : {fastq_dir}
    * GROUP FILE : {group_file}
    * REFERENCE FILE : {ref_file}
    * VECTOR FILE : {vec_file}
    * CPU : {cpu}
    
    """

        elif not os.path.exists(fastq_dir) or not os.path.isdir(fastq_dir):
            text = f"### WARNING : Directory doesn't exist or is not a directory : {fastq_dir} "

        elif not os.path.exists(group_file) or not os.path.isfile(group_file):
            text = f"### WARNING : File doesn't exist or is not a file : {group_file}"

        elif not os.path.exists(ref_file) or not os.path.isfile(ref_file):
            text = f"### WARNING : File doesn't exist : {ref_file} or is not a file"

        elif not os.path.exists(vec_file) or not os.path.isfile(vec_file):
            text = f"### WARNING : File doesn't exist : {vec_file} or is not a file"

        else:
            at = 'success'
            text = f"""
    ### All fields successfully filled 
    
    <hr>
    * PROJECT NAME : {project_name}
    * OUTPUT DIRECTORY : {output_dir}
    * FASTQ DIRECTORY : {fastq_dir}
    * GROUP FILE : {group_file}
    * REFERENCE FILE : {ref_file}
    * VECTOR FILE : {vec_file}
    * CPU : {cpu}
    
    """
            variables = {
                "project_name": project_name,
                "out_dir": output_dir,
                "ref_file": ref_file,
                "vec_file": vec_file,
                "group_file": group_file,
                "fastq_dir": fastq_dir,
                "cpu": cpu
            }
            # print(variables)

            file = open(file_load, 'wb')
            pickle.dump(variables, file)
            file.close()

        result.object = text.format(alert_type="success")
        return

    print_btn.param.watch(print_value, 'clicks')
    init_btn.param.watch(reinit_form, 'clicks')
    load_btn.param.watch(load_form, 'clicks')

    button = pn.Row(print_btn, load_btn, init_btn)
    row1 = pn.Row(project_namef, out_dirf)
    col1 = pn.Column(row1, fastq_dirf, group_filef, ref_filef, vec_filef, cpu_pf, button, result, width=800)

    # box
    tab = pn.WidgetBox('# INPUT FORM', col1, background='#E3ECF1')
    display(tab)

    return project_namef, out_dirf, ref_filef, vec_filef, group_filef, fastq_dirf, cpu_pf


def convert_mb(size):
    return round(size / 1000000, 2)


def dashboard_genome(reference_genome, output_dir):

    from Bio import SeqIO

    #  list to collect for each chromosome, chromosome size
    genome = []

    # read fasta genome reference
    FastaFile = open(reference_genome, 'r')
    for rec in SeqIO.parse(FastaFile, 'fasta'):
        seq_name = rec.id
        seq_len = len(rec.seq)
        # add chromosome name and size into genome list
        genome.append([seq_name, seq_len])
    FastaFile.close()

    # Create pandaframe from list
    df_genome = pd.DataFrame(genome, columns=["chr", "size", ])

    ### Generate ideogram file used to generate map of contigs anchored on genome
    df_ideogram = df_genome
    df_ideogram['start'] = 0
    df_ideogram['end'] = df_ideogram['size'] - 1
    df_ideogram['name'] = df_ideogram['chr']
    df_ideogram['gieStain'] = 'gvar'

    ### calculate genome size
    total_genome = df_genome['size'].sum()
    df_genome['size'] = round(df_genome['size'].div(1000000), 2)

    ### plot genome size of each chromosome
    # Create a visualization
    plt.figure(figsize=(8, 6))
    ax = sns.barplot(
        data=df_genome,
        x="chr", y="size"
    )
    ax.set(xlabel="Chr", ylabel="Size (Mb)")
    plt.close()
    stat = f"## Genome size (Mb) : {df_genome['size'].sum()}"

    # create a title for the dashboard
    dashboard_title = '# Reference genome file '

    dashboard = pn.Column(dashboard_title, stat, pn.Row(ax.figure, pn.panel(df_genome)),
                          sizing_mode='stretch_both', background='WhiteSmoke').servable()
    display(dashboard)

    # drop column size of df_ideogram and write the df into the ideogram file
    df_ideogram.drop('size', 1, inplace=True)
    df_ideogram.to_csv(os.path.join(output_dir, "ideogram.txt"), sep='\t', index=False, header=False)

    return total_genome


def dashboard_group(df):
    pn.extension()

    # create a class containing an species selector drop-down, various plots, and a data table output
    class GroupDashboard(param.Parameterized):
        # drop down selector widget containing the list of species
        # Animal = param.ObjectSelector(default='Goat', objects=list(df.Animal.unique()))
        Species = param.ObjectSelector(default=df.Species.unique().any(), objects=list(df.Species.unique()))

        # create data set containing only the data applicable to the species  in the drop down selector
        def get_data(self):
            class_df = df[(df.Species == self.Species)].copy()
            return class_df

        # seaborn box plot for the chosen animal
        def box_view(self):
            data = df
            plt.figure()  # figsize=(4, 3))
            ax = sns.countplot(y=data['Species'], palette="hls")
            plt.close()
            return ax.figure

        # table of data for the chosen animal
        def table_view(self):
            data = self.get_data()
            return data

    species_list = df.Species.unique().tolist()
    stat = f"## Samples : {df['sample'].count()} \n#### Group : {' '.join(species_list)}\n "

    # create an instance of the class
    rd = GroupDashboard(name='')

    # create a title for the dashboard
    dashboard_title = '# Population Group'

    # create some text describing the dashboard
    dashboard_desc = 'Some statistics and list sample by group '

    # create a dashboard, defining the layout as one column containing the
    # dashboard title, dashboard description, drop down selector,
    # box plot, and data table
    dashboard = pn.Column(dashboard_title,
                          # dashboard_desc,
                          pn.Row(stat, rd.box_view),  # box plot
                          rd.param,
                          pn.Row(rd.table_view, scroll=True, width=800, height=400),  # data table
                          sizing_mode='stretch_both', background='WhiteSmoke', scroll=True
                          )

    # show the dashboard with the data embedded ,
    # (for using in an html download of the notebook so that
    # no 'callback' is required from the browser to get the data)
    dashboard.embed()


def dashboard_fastq(csv, total_size, df_group):
    col_list = ["sample", "file", "reads", "len", "len_mean", "len_min", "qual_min", "qual_max", "qual_mean", "%A",
                "%C", "%G", "%T", "%N", "total_bases"]
    df_fastq_stat = pd.read_csv(csv, index_col=False, sep="\t", usecols=col_list)
    df_fastq_stat.sort_values(by=['sample', 'file'], inplace=True)

    df_stat_merged = pd.merge(left=df_group, right=df_fastq_stat, left_on='sample', right_on='sample')

    species_count = df_stat_merged.groupby('sample')['Species'].unique().value_counts()

    table_sp = f"""| Samples by group  |  |\n|:---|:---:|\n"""
    for key, item in species_count.items():
        table_sp += f"| {key[0]} : | {int(item)} samples |\n"

    sns.set_style("darkgrid")

    # Plot total bases pb by sample
    df_merge = df_stat_merged.groupby('sample')['total_bases'].sum()
    plt.figure()  # figsize=(10, 5))
    total = df_merge.plot(kind='bar', color='#66C2A5')
    total.set(xlabel="sample", ylabel="Total (Mb)")
    total.set_title("Total pb sequenced by sample")
    plt.close()

    total_desc = f"""### About population
* {len(df_stat_merged['sample'].unique())} samples
* Total base sequenced : {int(df_stat_merged['total_bases'].sum() / 1000000)} Mb"""

    df_grouped = df_stat_merged.groupby('Species')['total_bases'].sum().div(1000000)

    table = f"""| Pb sequenced by groups |  |\n|:---|:---:|\n"""
    for key, item in df_grouped.items():
        table += f"| {key} : | {int(item)} Mb |\n"

    # Coverage
    df_merge = df_stat_merged.groupby('sample')['total_bases'].sum().div(total_size)
    plt.figure()  # figsize=(10, 5))
    x = df_merge.plot(kind='bar', color='#66C2A5')
    x.set(xlabel="sample", ylabel="X")
    x.set_title("Sequencing coverage by sample")
    plt.close()

    # print(df_stat_group.groupby('sample')['Species'].unique().values)
    x_values = df_merge.values
    # df_Xgrouped = df_stat_group.groupby('Species')['total_bases'].sum().div(total_size)
    # print(df_Xgrouped)
    total_x = f"""### About Sequencing coverage
* x : {int(x_values.min())} - {int(x_values.max())} 
* mean x : {int(x_values.mean())} """

    # Plot length read
    # df_read = df_stat_merged[["file","len_min","len","len_mean"]]
    plt.figure()  # figsize=(10, 5))
    read = sns.scatterplot(x='file', y='len_mean', data=df_stat_merged, linewidth=2, palette="colorblind")
    read.set(xlabel="file", ylabel="read length (pb)")
    read.set_title("Read length over samples")
    plt.close()

    read_desc = f"""### About read length
* Length read: {df_stat_merged['len_min'].min()} - {df_stat_merged['len'].max()}
* Mean length read : {round(df_stat_merged['len_mean'].mean(), 2)}"""

    # Plot quality
    # df_qual = df_fstat_merged[["file","qual_min","qual_max","qual_mean"]]
    plt.figure()  # figsize=(10, 5))
    qual = sns.scatterplot(x='file', y='qual_mean', data=df_stat_merged, linewidth=2.5, palette="bright")
    qual.set(xlabel="file", ylabel="qual mean")
    qual.set_title("Quality mean over samples")
    plt.close()

    qual_desc = f"""### About sequencing quality
* Quality: {df_stat_merged['qual_min'].min()} - {df_stat_merged['qual_max'].max()}
* Mean quality : {round(df_stat_merged['qual_mean'].mean(), 2)}"""

    # Plot base ratio
    df_ratio = df_stat_merged[["file", "%A", "%C", "%G", "%T", "%N"]]
    plt.figure()  # figsize=(10, 5))
    base = sns.lineplot(x='file', y='value', hue='variable', data=pd.melt(df_ratio, 'file'),
                        linewidth=2, markersize=8, palette="hls", marker='o', linestyle='dashed')
    base.set(xlabel="file", ylabel="base ratio")
    base.set_title("Base ratio (A,T,C,G,N)")
    plt.legend('ACGTN', ncol=1, loc='upper left');
    plt.close()

    base_desc = "### About base ratio"

    df_short = df_stat_merged[
        ["sample", "file", "reads", "len_mean", "qual_mean", "%A", "%C", "%G", "%T", "%N", "total_bases"]]
    row1 = pn.Row()
    dashboard_title = '# Some statistics about fastq files'
    dashboard = pn.Column(dashboard_title,
                          pn.Row(total_desc, pn.Spacer(sizing_mode='stretch_both', width=100),
                                 pn.pane.Markdown(table_sp), width=800),
                          pn.Row(pn.pane.Markdown(table), pn.Spacer(sizing_mode='stretch_both', width=100),
                                 total.figure, width=800),
                          pn.Row(total_x, pn.Spacer(sizing_mode='stretch_both', width=100), x.figure, width=800),
                          pn.Row(read_desc, pn.Spacer(sizing_mode='stretch_both', width=100), read.figure, width=800),
                          pn.Row(qual_desc, pn.Spacer(sizing_mode='stretch_both', width=100), qual.figure, width=800),
                          pn.Row(base_desc, pn.Spacer(sizing_mode='stretch_both', width=100), base.figure, width=800),
                          pn.panel(df_short, width=800), sizing_mode='stretch_both', background='WhiteSmoke',
                          width=800).servable()

    display(dashboard)


def stat(df, col_stat):
    df_stat = df.agg({col_stat: ['mean', 'min', 'max']})
    return df_stat[col_stat]['min'], df_stat[col_stat]['max'], round(df_stat[col_stat]['mean'], 2)


def group_stat(df, col_group, col_stat):
    return df.groupby(col_group).agg({col_stat: ['mean', 'min', 'max']})


def format_stat(min, max, mean):
    return f"""<br />\n
* min - max : {min} - {max} %
* mean : {mean}"""



def format_stat_int(min, max, mean):
    return f"""<br />\n
* min - max : {int(min)} - {int(max)}
* mean : {int(mean)}"""



def dashboard_flagstat(stat_file, df_group):
    df_bam_stat = pd.read_csv(stat_file, index_col=False, sep=",")
    df_bam_stat.sort_values(by=['sample'], inplace=True)

    df_stat_merged = pd.merge(left=df_bam_stat, right=df_group, left_on='sample', right_on='sample')

    # check pop number
    species_count = df_stat_merged.groupby('sample')['Species'].unique().value_counts()
    table_sp = f"""| Samples by group  |  |\n|:---|:---:|\n"""
    for key, item in species_count.items():
        table_sp += f"| {key[0]} : | {int(item)} samples |\n"

    # Stat map
    mapped_min, mapped_max, mapped_mean = stat(df_stat_merged, 'MAPPED')
    mapped_title = "<br /><hr><br />\n### Reads mapped"
    mapped_desc = format_stat(mapped_min, mapped_max, mapped_mean)
    mapped_desc_sp = group_stat(df_stat_merged, 'Species', 'MAPPED')

    paired_min, paired_max, paired_mean = stat(df_stat_merged, 'PAIRED')
    paired_title = "<br /><hr><br />\n### Reads mapped properly paired"
    paired_desc = format_stat(paired_min, paired_max, paired_mean)
    paired_desc_sp = group_stat(df_stat_merged, 'Species', 'PAIRED')

    # stat unmapped
    unmapped_min, unmapped_max, unmapped_mean = stat(df_stat_merged, 'UNMAPPED')
    unmapped_title = "<br /><hr><br />\n### Reads unmapped"
    unmapped_desc = format_stat(unmapped_min, unmapped_max, unmapped_mean)
    unmapped_desc_sp = group_stat(df_stat_merged, 'Species', 'UNMAPPED')

    plt.figure()  # figsize=(12, 5))
    ratio = sns.scatterplot(x='sample', y='value', hue='variable', data=pd.melt(df_bam_stat, 'sample'))
    ratio.set(xlabel="sample", ylabel="ratio")
    ratio.set_title("Read mapped ratio ")
    plt.legend(bbox_to_anchor=(1.01, 1), borderaxespad=0);
    plt.close()

    from io import StringIO

    groups = pn.widgets.MultiChoice(
        name='Group', options=list(df_stat_merged.Species.unique()), margin=(0, 20, 0, 0)
    )
    samples = pn.widgets.MultiChoice(
        name='Sample', options=list(df_stat_merged['sample']), margin=(0, 20, 0, 0)
    )

    @pn.depends(groups, samples)
    def filtered_smp(grp, smp):
        df = df_stat_merged.sort_values(by=['sample'], inplace=True)
        if groups.value:
            df = df_stat_merged[df_stat_merged.Species.isin(grp)]
        if samples.value:
            df = df_stat_merged[df_stat_merged['sample'].isin(smp)]
        return df

    dashboard_title = '# Some statistics about mapping step'
    dashboard = pn.Column(dashboard_title,
                          pn.pane.Markdown(table_sp),
                          mapped_title,
                          pn.Row(mapped_desc, pn.panel(mapped_desc_sp, width=400)),
                          paired_title,
                          pn.Row(paired_desc, pn.panel(paired_desc_sp, width=400)),
                          unmapped_title,
                          pn.Row(unmapped_desc, pn.panel(unmapped_desc_sp, width=400)),
                          ratio.figure,
                          pn.Row(groups, samples),
                          pn.panel(filtered_smp, width=400, max_height=400, scroll=True),
                          sizing_mode='stretch_both', background='WhiteSmoke', width=800).servable()
    display(dashboard)


def box_config_abyss(df):
    # cmd result
    text = "No filled"
    at = 'warning'
    result = pn.pane.Alert(text.format(alert_type=at), alert_type=at, height=200)

    df_sorted = df.sort_values(by=['sample'], ascending=True)

    # buton
    print_btn = pn.widgets.Button(name='SAVE', width=100, button_type='primary')
    init_btn = pn.widgets.Button(name='INIT', width=100, button_type='primary')

    # form
    k = pn.widgets.RangeSlider(name='K-mer length', start=24, end=92, step=4, value=(64, 68))
    step = pn.widgets.IntInput(name='Step', value=4, step=1, start=1, end=10)
    threshold = pn.widgets.IntInput(name='Minimal length to filter', value=300, step=100, start=100, end=10000)

    accession = pn.widgets.MultiSelect(name='Accession', options=list(df_sorted['sample']),
                                       value=[df_sorted['sample'][0], df_sorted['sample'][1]], size=4)

    def print_value(event):
        at = 'success'
        text = f"""
            ### k successfully filled

            <hr>
             * k : {k.value}
             * step : {step.value}
             * Accession : {accession.value}
             * Minimal length : {threshold.value}
        """

        result.object = text.format(alert_type="success")
        return

    def reinit_form(event):
        k.value = (64, 68)
        accession.value = list()
        step.value = 4
        threshold.value = 300
        text = "empty"
        return

    print_btn.param.watch(print_value, 'clicks')
    init_btn.param.watch(reinit_form, 'clicks')

    button = pn.Row(print_btn, init_btn)
    row1 = pn.Row(k, step)
    col1 = pn.Column(row1, accession, threshold, button, result, width=800)

    # box
    tab = pn.WidgetBox('# INPUT FORM', col1, background='#E3ECF1')
    display(tab)

    return k, step, accession, threshold


def box_config_abyss2():
    # cmd result
    text = "No filled"
    at = 'warning'
    result = pn.pane.Alert(text.format(alert_type=at), alert_type=at, height=200)

    # buton
    print_btn = pn.widgets.Button(name='SAVE', width=100, button_type='primary')
    init_btn = pn.widgets.Button(name='INIT', width=100, button_type='primary')

    # form
    k = pn.widgets.IntInput(name='K-mer length', start=24, end=92, step=1, value=64)

    def print_value(event):
        at = 'success'
        text = f"""
            ### k successfully filled

            <hr>
             * k : {k.value}
        """

        result.object = text.format(alert_type="success")
        return

    def reinit_form(event):
        k.value = 64
        text = "empty"
        return

    print_btn.param.watch(print_value, 'clicks')
    init_btn.param.watch(reinit_form, 'clicks')

    button = pn.Row(print_btn, init_btn)
    col1 = pn.Column(k, button, result, width=800)

    # box
    tab = pn.WidgetBox('# INPUT FORM', col1, background='#E3ECF1')
    display(tab)

    return k


def dashboard_ab(stat_l, stat_N, stat_L, diri):
    stat_len_df = pd.read_csv(os.path.join(diri, "assembly-stats-" + stat_l[0] + ".csv"), sep='\t')

    sns.set_style("darkgrid")

    table_ln = f"""| Total length assembled  |  |\n|:---|:---:|\n"""
    plt.figure()  #
    size_ln = sns.relplot(x='k', y='value', hue='stat', data=stat_len_df, col='id', kind="line")
    size_ln.set(ylabel="Ln", xlabel="kmer")
    plt.close()

    table_n = f"""| N assembled  |  |\n|:---|:---:|\n"""
    stats_N_files = []
    for stat in stat_N:
        stats_N_files.append(os.path.join(diri, "assembly-stats-" + stat + ".csv"))

    stats_N_df = pd.concat([pd.read_csv(f, sep='\t') for f in stats_N_files], ignore_index=True)
    plt.figure()
    size_N = sns.relplot(x='k', y='value', hue='stat', data=stats_N_df, col='id', kind="line")
    size_N.set(ylabel="N", xlabel="kmer")
    plt.close()

    table_l = f"""| TL assembled  |  |\n|:---|:---:|\n"""
    stats_L_files = []
    for stat in stat_L:
        stats_L_files.append(os.path.join(diri, "assembly-stats-" + stat + ".csv"))
    stats_L_df = pd.concat([pd.read_csv(f, sep='\t') for f in stats_L_files], ignore_index=True)

    plt.figure()
    size_L = sns.relplot(x='k', y='value', hue='stat', data=stats_L_df, col='id', kind="line")
    size_L.set(ylabel="L", xlabel="kmer")
    plt.close()
    # size_L.savefig(diri + "stat_L.png")

    dashboard_title = '# Some statistics about abyss assembly'
    dashboard = pn.Column(dashboard_title,
                          pn.Row(pn.pane.Markdown(table_ln)),
                          pn.Row(size_ln.figure, width=800),
                          pn.Row(pn.pane.Markdown(table_n)),
                          pn.Row(size_N.figure, width=800),
                          pn.Row(pn.pane.Markdown(table_l)),
                          pn.Row(size_L.figure, width=800),
                          sizing_mode='stretch_both', background='WhiteSmoke', width=800).servable()
    display(dashboard)


def dashboard_assembly(stat_file, df_group):

    # read stat file
    df_stat = pd.read_csv(stat_file, index_col=False, sep="\t")

    # extract  sample name from the column "filename"
    df_stat['sample'] = df_stat.filename.apply(lambda x: os.path.basename(x))
    df_stat['sample'] = df_stat['sample'].str.split('_', expand=True)[0]

    # Merge the dataframes "stat" and "group" into a new dataframe
    df_stat_merged = pd.merge(left=df_stat, right=df_group, left_on='sample', right_on='sample')
    df_stat_merged = df_stat_merged[
        ["sample", "number", "mean_length", "shortest", "longest", "total_length", "Species"]]
    df_stat_merged.sort_values(by=['sample'], inplace=True)

    ctg_min, ctg_max, ctg_mean = stat(df_stat_merged, 'number')
    pb_min, pb_max, pb_mean = stat(df_stat_merged, 'total_length')
    ctg_pb_mean = df_stat_merged.mean_length.mean()
    ctg_pb_max = df_stat_merged.longest.max()
    ctg_pb_min = df_stat_merged.shortest.min()
    df_stat_merged = df_stat_merged.rename(
        columns={"number": "Contigs number", "shortest": "Shortest (pb)", "longest": "Longuest (pb)",
                 "total_length": "total_length (pb)"})

    # Stat contigs numbers
    ctg_title = "<br /><hr><br />\n### Contigs number"
    ctg_desc = format_stat_int(ctg_min, ctg_max, ctg_mean)
    ctg_desc_sp = group_stat(df_stat_merged, 'Species', 'Contigs number')

    # Stat contigs length
    ctg_pb_title = f"<br /><hr><br />\n### Contigs size"
    ctg_pb_desc = format_stat(convert_mb(pb_min), convert_mb(pb_max), convert_mb(pb_mean))
    ctg_pb_sp = group_stat(df_stat_merged, 'Species', 'total_length (pb)')

    # Stat contogh length  ctg_pb_mean
    ctg_ln_title = f"<br /><hr><br />\n### Contigs length"
    ctg_ln_desc = format_stat_int(ctg_pb_min, ctg_pb_max, ctg_pb_mean)

    # plt.figure()#figsize=(12, 5))
    # ratio = sns.scatterplot(x='sample', y='value', hue='variable', data=pd.melt(df_stat_merged, 'sample'))
    # ratio.set(xlabel="sample", ylabel="Contigs number")
    # ratio.set_title("Contigs number")
    # plt.legend(bbox_to_anchor=(1.01, 1), borderaxespad=0);
    # plt.close()

    from io import StringIO

    groups = pn.widgets.MultiChoice(
        name='Group', options=list(df_stat_merged.Species.unique()), margin=(0, 20, 0, 0)
    )
    samples = pn.widgets.MultiChoice(
        name='Sample', options=list(df_stat_merged['sample']), margin=(0, 20, 0, 0)
    )

    @pn.depends(groups, samples)
    def filtered_smp(grp, smp):
        df = df_stat_merged.sort_values(by=['sample'], inplace=True)
        if groups.value:
            df = df_stat_merged[df_stat_merged.Species.isin(grp)]
        if samples.value:
            df = df_stat_merged[df_stat_merged['sample'].isin(smp)]
        return df

    dashboard_title = '# Some statistics about abyss step'
    dashboard = pn.Column(dashboard_title,
                          ctg_title,
                          pn.Row(ctg_desc, pn.panel(ctg_desc_sp, width=400)),
                          ctg_pb_title,
                          pn.Row(ctg_pb_desc, pn.panel(ctg_pb_sp, width=400)),
                          ctg_ln_title,
                          pn.Row(ctg_ln_desc),
                          pn.Row(groups, samples),
                          pn.panel(filtered_smp, width=800, max_height=600, scroll=True),
                          sizing_mode='stretch_both', background='WhiteSmoke', width=800).servable()
    display(dashboard)


def dashboard_ass(stat_file, df_group):
    # read stat file
    df_stat = pd.read_csv(stat_file, index_col=False, sep="\t")
    # print(df_stat)

    ctg_num = df_stat.number.iloc[0]
    pb_tot = df_stat.total_length.iloc[0]
    ctg_pb_mean = df_stat.mean_length.iloc[0]
    ctg_pb_max = df_stat.longest.iloc[0]
    ctg_pb_min = df_stat.shortest.iloc[0]

    #print(ctg_pb_max)
    # Stat contigs numbers
    ctg_title = f"<br /><hr><br />\n### Contigs number : {ctg_num} ({convert_mb(pb_tot)} Mb)"

    # Stat contogh length  ctg_pb_mean
    ctg_ln_title = f"<br /><hr><br />\n### Contigs length : "
    ctg_ln_desc = format_stat_int(ctg_pb_min, ctg_pb_max, ctg_pb_mean)

    dashboard_title = '# Some statistics'
    dashboard = pn.Column(dashboard_title,
                          ctg_title,
                          pn.Row(ctg_ln_desc),
                          sizing_mode='stretch_both', background='WhiteSmoke', width=800).servable()
    display(dashboard)


def dashboard_cdhit(df_cdhit):
    df_cdhit['pb'] = df_cdhit['pb'].astype(int)

    cdhit_comp = f"""### Cluster composition
* Sequences Number : {df_cdhit['pb'].count()} ({df_cdhit['pb'].sum()} pb)
    * singleton: {df_cdhit[df_cdhit.ln == 1]['pb'].count()}
    * clusters: {df_cdhit[(df_cdhit.ln > 1)]['index'].count()}
"""

    cdhit_desc = f"""### {df_cdhit['pb'].count()} sequences 
* min -  max length : {df_cdhit['pb'].min()} bp - {df_cdhit['pb'].max()} bp
* mean length : {int(df_cdhit['pb'].mean())} bp
"""

    singleton_desc = f"""### {df_cdhit[df_cdhit.ln == 1]['pb'].count()} singletons 
* min - max length : {df_cdhit[(df_cdhit.ln == 1)].pb.min()} bp -  {df_cdhit[(df_cdhit.ln == 1)].pb.max()} bp
* mean length : {int(df_cdhit[(df_cdhit.ln == 1)].pb.mean())} bp
"""
    cluster_desc = f"""### {df_cdhit[(df_cdhit.ln > 1)]['index'].count()} clusters
* min - max length : {df_cdhit[(df_cdhit.ln > 1)].pb.min()} bp - {df_cdhit[(df_cdhit.ln > 1)].pb.max()} bp
* mean length : {int(df_cdhit[(df_cdhit.ln > 1)].pb.mean())}
* from {df_cdhit[(df_cdhit.ln > 1)].ln.min()} to  {df_cdhit[(df_cdhit.ln > 1)].ln.max()}
* with ~ {int(df_cdhit[(df_cdhit.ln > 1)].ln.mean())} sequences
"""

    # plt.figure(figsize=(17, 6))
    plt.figure(figsize=(10, 5))
    nb = sns.countplot(y=df_cdhit[(df_cdhit.ln < 20)].ln, palette="hls")
    nb.set(xlabel="cluster number", ylabel="Sequence number within a cluster")
    nb.set_title("Distribution of sequence number by cluster (< 20 sequences within a cluster)")
    plt.close()

    plt.figure(figsize=(10, 5))
    size = sns.histplot(data=df_cdhit[(df_cdhit.pb < 3000)], x='pb')
    size.set(ylabel="Sequence number", xlabel='Length')
    size.set_title("Distribution of sequence length")
    plt.close()

    # plt.figure(figsize=(17, 6))
    # ax = plt.gca()
    plt.figure(figsize=(10, 5))
    sc = sns.scatterplot(x='pb', y='ln', data=df_cdhit, palette='blue')
    sc.set(xlabel="pb", ylabel="Cluster sequence number")
    sc.set_title("ength of the contigs versus nomber of sequences within a cluster")
    plt.close()

    dashboard_title = '# Some statistics about clustering step'
    dashboard = pn.Column(dashboard_title,
                          pn.pane.Markdown(cdhit_desc),
                          pn.Row(pn.pane.Markdown(cluster_desc), pn.pane.Markdown(singleton_desc)),
                          nb.figure,
                          size.figure,
                          sc.figure,
                          # pn.Row(sin_ln.figure, cl_ln.figure),
                          # pn.panµel(filtered_smp, width=400, max_height=400, scroll=True),
                          sizing_mode='stretch_both', background='WhiteSmoke', width=800).servable()
    display(dashboard)




# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for genes).  The rest of this script will be prepping data
# for input to this function
# Source : PUT SOURCE
def chromosome_collections(df, y_positions, height, **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        # print(chrom)
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


def plot_map(output_dir, bed_file):
    ###### Fig parameters
    chrom_height = 0.5  # Height of each ideogram
    chrom_spacing = 0.8  # Spacing between consecutive ideograms
    gene_height = 0.3  # Height of the gene track. Should be smaller than `chrom_spacing` in order to fit correctly
    gene_padding = 0.1  # Padding between the top of a gene track and its corresponding ideogram

    figsize = (14, 16)  # Width, height (in inches)

    #####  Read ideogram.txt
    ideo = pd.read_table(os.path.join(output_dir, "ideogram.txt"),
                             names=['chrom', 'start', 'end', 'name', 'gieStain']
                             )

    # Add a new column for width
    ideo['width'] = ideo.end - ideo.start

    chromosome_list = ideo['chrom'].to_list()

    # Keep track of the y positions for ideograms and genes for each chromosome,
    # and the center of each ideogram (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    gene_ybase = {}
    chrom_centers = {}

    # Iterate in reverse so that items in the beginning of `chromosome_list` will
    # appear at the top of the plot
    for chrom in chromosome_list[::-1]:
        # print(chrom)
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        gene_ybase[chrom] = ybase - gene_height - gene_padding
        ybase += chrom_height + chrom_spacing

    # Colors for different chromosome stains
    color_lookup = {
        'gneg': (1., 1., 1.),
        'gpos25': (.6, .6, .6),
        'gpos50': (.4, .4, .4),
        'gpos75': (.2, .2, .2),
        'gpos100': (0., 0., 0.),
        'acen': 'grey',  # "steelblue", #(.8, .4, .4),
        'gvar': '#f7dd8d',  # "steelblue", #(.8, .8, .8),
        'stalk': '#f7dd8d'  # "lightblue", #(.9, .9, .9)}
    }

    # Add a new column for colors
    ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])

    # Same thing for genes
    genes = pd.read_table(
        bed_file, header=1,
        names=['chrom', 'CTG_name', 'start', 'end'],
        usecols=range(4))

    genes["start"] = pd.to_numeric(genes["start"])
    genes["end"] = pd.to_numeric(genes["end"])
    genes = genes[genes.chrom.apply(lambda x: x in chromosome_list)]
    genes['width'] = genes.end - genes.start
    genes['colors'] = 'red'  # 2243a8'

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # Now all we have to do is call our function for the ideogram data...
    # print("adding ideograms...")
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height):
        # print(collection)
        ax.add_collection(collection)

    # ...and the gene data
    # print("adding genes...")
    for collection in chromosome_collections(genes, gene_ybase, gene_height, alpha=0.5, linewidths=0):
        ax.add_collection(collection)

    # Axes tweaking
    ax.set_facecolor('white')  ##fff2cc')
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list, fontsize=12)
    ax.axis('tight')

    panref_dir = os.path.dirname(bed_file)
    plt.savefig(os.path.join(panref_dir, 'contigs_on_genome_map.png'))
    plt.close(fig)


def dashboard_anchoring(ctg_fasta, panref_keep_file, panref_bed_file, output_dir, anc_stat_dict):
    # Contigs number total
    nb_tot_ctgs = int(get_seq_nb(ctg_fasta))

    # stats
    # print(anc_stat_dict)
    df_thebest = pd.read_csv(panref_keep_file, sep=",")

    # print(df_thebest)
    anc_comp = f"""### Anchoring process
* {anc_stat_dict['contigs_anchored']} contigs anchored over the {nb_tot_ctgs} contigs analyzed 
* {anc_stat_dict['uniq_anchored']} contigs placed at an unique position  with a sequence length mean of  {anc_stat_dict['uniq_anchored_mean']} pb
    """

    anc_desc = f"""### {anc_stat_dict['uniq_pass']} contigs pass the filter of length pb  - {round(anc_stat_dict['uniq_pass'] / nb_tot_ctgs * 100, 2)} % of all ctgs
* min -  max length : {df_thebest.CTG_len.min()} bp - {df_thebest.CTG_len.max()}  bp
* mean length : {int(df_thebest.CTG_len.mean())}  bp
"""
    # print(df_thebest.sort_values("CHR_name").head())
    df_thebest = df_thebest.sort_values("CHR_name")
    # print(df_thebest.head())

    table_sp = f"""| Chromosome | Contigs anchored |\n|:---|:---:|\n"""
    for key, item in df_thebest.CHR_name.value_counts().items():
        table_sp += f"""| {key}  | {int(item)} |\n"""

    #     singleton_desc = f"""### {df_cdhit[df_cdhit.ln == 1]['pb'].count()} singletons
    # * min - max length : {df_cdhit[(df_cdhit.ln == 1)].pb.min()} bp -  {df_cdhit[(df_cdhit.ln == 1)].pb.max()} bp
    # * mean length : {int(df_cdhit[(df_cdhit.ln == 1)].pb.mean())} bp
    # """
    #     cluster_desc = f"""### {df_cdhit[(df_cdhit.ln > 1)]['index'].count()} clusters
    # * min - max length : {df_cdhit[(df_cdhit.ln > 1)].pb.min()} bp - {df_cdhit[(df_cdhit.ln > 1)].pb.max()} bp
    # * mean length : {int(df_cdhit[(df_cdhit.ln > 1)].pb.mean())}
    # * from {df_cdhit[(df_cdhit.ln > 1)].ln.min()} to  {df_cdhit[(df_cdhit.ln > 1)].ln.max()}
    # * with ~ {int(df_cdhit[(df_cdhit.ln > 1)].ln.mean())} sequences
    # """

    # plt.figure(figsize=(17, 6))
    plt.figure(figsize=(10, 5))
    # sns.countplot(x ='sex', data = df)
    nb = sns.countplot(x=df_thebest.CHR_name, palette="hls")
    nb.set(xlabel="chromosome", ylabel="Number of contigs anchored")
    nb.set_title("Distribution of contigs by chromosome ")
    plt.close()

    #     plt.figure(figsize=(10, 5))
    #     size = sns.histplot(data=df_cdhit[(df_cdhit.pb < 3000)], x='pb')
    #     size.set(ylabel="Sequence number", xlabel='Length')
    #     size.set_title("Distribution of sequence length")
    #     plt.close()

    #     # plt.figure(figsize=(17, 6))
    #     # ax = plt.gca()
    #     plt.figure(figsize=(10, 5))
    #     sc = sns.scatterplot(x='pb', y='ln', data=df_cdhit, palette='blue')
    #     sc.set(xlabel="pb", ylabel="Cluster sequence number")
    #     sc.set_title("ength of the contigs versus nomber of sequences within a cluster")
    #     plt.close()

    plot_map(output_dir, panref_bed_file)
    panref_dir = os.path.dirname(panref_bed_file)
    panref_map_file = os.path.join(panref_dir, 'contigs_on_genome_map.png')
    dashboard_title = '# Some statistics about anchoring step'
    dashboard = pn.Column(dashboard_title,
                          pn.pane.Markdown(anc_comp),
                          pn.Row(pn.pane.Markdown(anc_desc, width=1000)),
                          pn.Row(nb.figure, pn.pane.Markdown(table_sp, width=200)),
                          pn.pane.PNG(panref_map_file),
                          # size.figure,
                          # sc.figure,
                          # pn.Row(sin_ln.figure, cl_ln.figure),
                          # pn.panµel(filtered_smp, width=400, max_height=400, scroll=True),
                          sizing_mode='stretch_both', background='WhiteSmoke', width=1500).servable()
    display(dashboard)

def get_seq_nb(fasta):

    result = subprocess.run(['grep', '-c', '>', fasta], capture_output=True, text=True)
    return (result.stdout)

# Dashboard