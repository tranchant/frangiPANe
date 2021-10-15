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
from ipywidgets import Button, Layout
import traitlets
from IPython.display import display,HTML
from tkinter import Tk, filedialog
#from IPython.core.display import HTML

import pandas as pd
import logging,os
from time import localtime, strftime


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
    return logger
    
def hide_cell():

    tag = HTML('''<script>
    code_show=true; 
    function code_toggle() {
        if (code_show){
            $('div.cell.code_cell.rendered.selected div.input').hide();
        } else {
            $('div.cell.code_cell.rendered.selected div.input').show();
        }
        code_show = !code_show
    } 

    $( document ).ready(code_toggle);
    </script>

    <a href="javascript:code_toggle()">[raw code]</a>''')
    display(tag)

def add_css() : 
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

def msg_button(desc,button_col, text_col):

    color_text={
        'warning' : 'mywarning',
        'classic' : 'myclassic'
    }

    color_button={
        'warning' : '#EC4646',
        'yellow' : '#F8DC81',
        'green': '#9FD89F',
        'blue' : '#9FD8DF'
    }


    btn_one = widgets.Button(description=desc, style={'button_color': color_button[button_col]}, layout=Layout(width='100%', height='40px'))
    btn_one.add_class(color_text[text_col])

    return(widgets.VBox([btn_one]))

def box_config():
#CD : PUT URL CODE ORIGINE
    class SelectDirButton(widgets.Button):
        """A file widget that leverages tkinter.filedialog."""

        def __init__(self) -> object:
            super(SelectDirButton, self).__init__()
            # Add the selected_files trait
            self.add_traits(dirs=traitlets.traitlets.List())
            # Create the button.
            self.description = "Select dir"
            self.icon = "square-o"
            self.style.button_color = "#E7E3D5"
            self.layout=Layout(width='40%', height='40px')
            # Set on click behavior.
            self.on_click(self.select_dir)

        @staticmethod
        def select_dir(b):
            """Generate instance of tkinter.filedialog.

            Parameters
            ----------
            b : obj:
                An instance of ipywidgets.widgets.Button
            """
            with outd:
                try:
                    # Create Tk root
                    root = Tk()
                    # Hide the main window
                    root.withdraw()
                    # Raise the root to the top of all windows.
                    root.call('wm', 'attributes', '.', '-topmost', True)
                    # List of selected fileswill be set to b.value
                    b.dir = filedialog.askdirectory()

                    b.description = os.path.basename(b.dir)
                    b.icon = "check-square-o"
                    b.style.button_color = "#9FD89F"
                except:
                    pass

    class SelectFilesButton(widgets.Button):
        """A file widget that leverages tkinter.filedialog."""

        def __init__(self) -> object:
            super(SelectFilesButton, self).__init__()
            # Add the selected_files trait
            self.add_traits(files=traitlets.traitlets.List())
            # Create the button.
            self.description = "Select file"
            self.icon = "square-o"
            self.style.button_color = "#E7E3D5"
            self.layout=Layout(width='20%', height='40px')
            # Set on click behavior.
            self.on_click(self.select_file)

        @staticmethod
        def select_file(b):
            """Generate instance of tkinter.filedialog.

            Parameters
            ----------
            b : obj:
                An instance of ipywidgets.widgets.Button
            """
            with out:
                try:
                    # Create Tk root
                    root = Tk()
                    # Hide the main window
                    root.withdraw()
                    # Raise the root to the top of all windows.
                    root.call('wm', 'attributes', '.', '-topmost', True)
                    # List of selected fileswill be set to b.value
                    b.file = filedialog.askopenfilename(multiple=False)

                    b.description = os.path.basename(b.file)
                    b.icon = "check-square-o"
                    b.style.button_color = "#9FD89F"
                except:
                    pass

    html1_value = widgets.HTML(
        value="<b><center><FONT face = 'WildWest' color='#255E7A'>INPUT & OUTPUT DATA</FONT></center></b>",
        placeholder='Some HTML',
        description='',
    )
    ##### DATA CONFIG
    style = {'description_width': 'initial'}
    project = widgets.Text(
        value='test',
        description='Project name:',
        layout=Layout(width='40%', height='40px'),
        style=style
    )

    outd = widgets.Output()
    label_dirw = widgets.Label("Output directory")
    out_dirw = SelectDirButton()
    out_dirw.description = "Output dir"

    out = widgets.Output()
    label_ref = widgets.Label("Genome Reference")
    ref_fastaw = SelectFilesButton()
    ref_fastaw.description = "Genome reference"

    label_fastq = widgets.Label("Fastq directory")
    fastq_dirw = SelectDirButton()
    fastq_dirw.description = "Fastq directory"

    label_group = widgets.Label("Group file")
    group_filew = SelectFilesButton()
    group_filew.description = "Group file"

    #index = widgets.Checkbox(
    #       description='Indexation?')

    top_box = widgets.HBox([label_dirw, outd, out_dirw, label_group, out, group_filew])
    bottom_box = widgets.HBox([label_fastq,outd, fastq_dirw, label_ref, out, ref_fastaw ])
    #box_data = widgets.VBox([html1_value, project, outd, out_dirw, out, ref_fastaw, group_filew, outd, fastq_dirw])
    box_data = widgets.VBox([html1_value, project, top_box, bottom_box])

    ##### CONFIG BOX
    df = pd.read_csv("conf/cluster.conf",names=["tag","val"], sep="=")
    # widget_display(df)
    cpu_max = int(df.iat[0, 1])
    partition = df.iat[1, 1].split(",")

    cpu = widgets.BoundedIntText(
        value=6,
        min=1,
        max=cpu_max,
        step=1,
        description='CPUs:',
        disabled=False
    )

    queue = widgets.Select(
        options=partition,
        description='Slurm partition:',
        disabled=False
    )

    slurm = widgets.Checkbox(
           description='Slurm',)

    box_cluster = widgets.VBox([cpu, queue, slurm])

    #### MAPPING CONFIG
    box_mapping = widgets.VBox([cpu])


    #### VBOX
    # defining a list with the contents of our windows
    children = [box_data, box_mapping] # box_cluster,
    # initializing a tab
    tab = widgets.Tab()
    # setting the tab windows
    tab.children = children
    # changing the title of the first and second window
    tab.set_title(0, 'Data')
    #tab.set_title(1, 'Cluster')
    tab.set_title(1, 'CPU')
    display(tab)
    return project, out_dirw, ref_fastaw, group_filew, fastq_dirw, cpu #, queue, slurm  #out_dirw,, fastq_dirw, group_filew,

