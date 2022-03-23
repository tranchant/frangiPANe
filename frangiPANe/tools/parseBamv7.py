#!/usr/bin/env python3

# License GNU/GPL
# Intellectual property belongs to IRD
# written by Christine Tranchant-Dubreuil (UMR DIADE/itrop platform - IRD)

"""

    Use it to place contigs on chromosomes using reads mapping information
    :Example:
    >>> python ~/Script/Python/parseBamv7.py -b bam_dir -d 10 -o file

"""

import argparse
import os.path
import pysam

from numpy import nanmean, nanstd, mean, std, abs, percentile


### IQR tells how spread the middle values are. It can be used to tell when a value is too far
### from the middle.
### An outlier is a point which falls more than 1.5 times the interquartile range above the third
### quartile or below the first quartile.
def iqr(dataset):

    sorted(dataset)

    # third and first quartile
    q1, q3 = percentile(dataset, [25, 75])

    # Find IQR which is the difference between third and first quartile
    iqr = q3 - q1

    # Find lower and upper bound
    lower_bound = q1 - (1.5 * iqr)
    upper_bound = q3 + (1.5 * iqr)

    return(lower_bound,upper_bound)




### read a list of bam file
### Store in a dictionary all informations related to the ctgs placement if
###     - one read maps in the 300 first or end pb
###     - the mate read maps on a chromosome
def read_bam(bam_list, depth, out_file):

    # chrList =  ['Chr01','Chr02','Chr03','Chr04','Chr05','Chr06','Chr07','Chr08','Chr09','Chr10','Chr11','Chr12','ChrUN-Ctg37','ChrUN-Ctg38','ChrUN-Ctg39','ChrUN-Ctg40','ChrUN-Ctg41','ChrUN-Ctg42','ChrUN-Ctg43','ChrUN-Ctg44','ChrUN-Ctg45','ChrUN-Ctg46','ChrUN-Ctg47','ChrUN-Ctg48','ChrUN-Ctg49','ChrUN-Ctg50','ChrUN-Ctg51','ChrUN-Ctg52','ChrUN-Ctg53','ChrUN-Ctg54']

    ### dictionary initialisation
    myRef = { }         # store informations about the chromosome and their length
    anchor = { }        # store all informations about the ctgs placement on chromosomes
    # Anchor has 2 firsts keys, 1rst the contig name and 2nd the chromosome name
    # And after a dictionay with all informations about the placement after reading successivly of all bam files
    # anchor[ctg_name][chr_name] =
    # {
         # 'posi': [(read_start on chromosome, mate_read_start on contig)],  # chr and ctg positions list for each read pairs
         # 'posiChr': [read_start],       # list of positions on chromosome
         # 'posiCtg': [read2_start],      # list of positions on contigs
         # 'rgList': [],                  # list of samples used to the placement
         # 'totalHit': 0,                 # Total number of hit along the contig
         # 'startHit': 0,                 # Total number of hit in the 300 first base pairs of the ctg
         # 'startPosiChr': [],            # Chr positions list when the mate read maps in the 300st pb of the ctg
         # 'startPosiCtg': [],            # Ctg positions list when the read maps in the 300st pb of the ctg
         # 'endHit': 0,                   # Total number of hit in the 300 last base pairs of the ctg
         # 'endPosiChr': [],              # Chr positions list when the mate read maps in the 300 last pb of the ctg
         # 'endPosiCtg': [],              # Ctg positions list when the read maps in the 300 last pb of the ctg
         # 'len': myRef[read_ref_name]    # Chromosome length
    # }

    #######################################"
    # Parsing bam files given in argument
    for bam_file in bam_list:

        if not '.bam' in bam_file:
            continue

        # get the sample name (TOGGLE nomenclature)
        # PATH/bam_dir/AAOSW.PICARDTOOLSSORT.bam => rg = AAOWS
        rg = bam_file.split("/")[-1].split(".")[0]

        # next file if bai file
        if '.bai' in bam_file:
            continue

        # use pysam to read bam file
        print(bam_file)
        print(f"--{rg}--")
        bamfile = pysam.AlignmentFile(bam_file, "rb")
        print(bamfile.mapped, bamfile.unmapped)

        # Extract reference name and length from the header
        if (not myRef):
            myRef = dict([(ref['SN'], ref['LN']) for ref in bamfile.header['SQ']])

        # variable inisialisation
        nb = 0              # total reads number in the bam file
        nbNotProp = 0       # count the number of reads that are not propoerly mapped
        nbNotProp2 = 0      # count the number of reads (and mate) placed on chr (and ctg) _ CASE1
        nbNotProp3 = 0      # count the number of reads (and mate) placed on chr (and ctg) _ CASE2
        nbUnmapped = 0
        nbMateUnmapped = 0
        nbSup = 0           # number of secondary reads
        nbSupTot = 0        # number of supplementary reads
        nb_noqual = 0

        # read bamfile line by line
        for read in bamfile:
            qual = 1
            if read.mapping_quality < 10:
                #print(read.reference_name, read.mapping_quality )
                nb_noqual += 1
                qual = 0
                continue

            # just read line if read is not proper mapped and is read forward
            if not read.is_proper_pair and read.is_read1:

                nbNotProp += 1

                # Exclude read and mate read not mapped
                if read.reference_name is None and read.next_reference_name is None:
                    nbNotProp -= 1
                    continue

                # Get the reference on which the read is mapped (read1)
                # and the reference on which its mate is mapped (read2)
                read1_ref_name = read.reference_name
                read2_ref_name = read.next_reference_name

                # CASE 1: read1 maps on a chromosome and mate read on a contig
                if 'Chr' in read1_ref_name and not 'Chr' in read2_ref_name:

                    # Get the start position of read1 on chromosome
                    read1_posi_list = read.get_reference_positions(full_length=True)
                    read1_start = read1_posi_list[0]

                    # Get the start position of mate read on the contig
                    read2_start = read.next_reference_start

                    # First time that this contig is found
                    if read2_ref_name not in anchor:
                        anchor[read2_ref_name] = {}

                    # First time that the chromosome is found with the contig
                    if read1_ref_name not in anchor[read2_ref_name]:
                        anchor[read2_ref_name][read1_ref_name] = {
                                                    'posi' :[(read1_start,read2_start)],
                                                    'posiChr': [read1_start],
                                                    'posiCtg': [read2_start],
                                                    'rgList' : [ ],
                                                    'totalHit': 0,
                                                    'startHit': 0,
                                                    'startPosiChr' : [],
                                                    'startPosiCtg': [],
                                                    'endHit' : 0,
                                                    'endPosiChr' : [ ],
                                                    'endPosiCtg': [],
                                                    'len' : myRef[read2_ref_name],
                                                    'qual' : [ ]
                                                    }
                    # The ctg placement on this chr has already been found precedently
                    # add the reads positions list on the ctg and chr in the dictionary
                    else:
                        anchor[read2_ref_name][read1_ref_name]['posi'].append((read1_start,read2_start))
                        anchor[read2_ref_name][read1_ref_name]['posiChr'].append(read1_start)
                        anchor[read2_ref_name][read1_ref_name]['posiCtg'].append(read2_start)
                        anchor[read2_ref_name][read1_ref_name]['qual'].append(qual)

                    nbNotProp2 += 1
                    anchor[read2_ref_name][read1_ref_name]['totalHit'] += 1

                    # If read2 mapped in the 300 st pb of the contigs
                    # add ctg and chr positions in "start" keys and the sample name
                    if read2_start is not None:
                        if (read2_start < 300):
                            anchor[read2_ref_name][read1_ref_name]['startHit'] += 1
                            anchor[read2_ref_name][read1_ref_name]['startPosiChr'].append(read1_start)
                            anchor[read2_ref_name][read1_ref_name]['startPosiCtg'].append(read2_start)
                            if rg not in anchor[read2_ref_name][read1_ref_name]['rgList']:
                                anchor[read2_ref_name][read1_ref_name]['rgList'].append(rg)

                        # If read2 mapped in the 300 last pb of the contigs
                        # add ctg and chr positions in "end" keys and the sample name
                        elif (read2_start > myRef[read2_ref_name] - 300):
                            anchor[read2_ref_name][read1_ref_name]['endHit'] += 1
                            anchor[read2_ref_name][read1_ref_name]['endPosiChr'].append(read1_start)
                            anchor[read2_ref_name][read1_ref_name]['endPosiCtg'].append(read2_start)
                            if rg not in anchor[read2_ref_name][read1_ref_name]['rgList']:
                                anchor[read2_ref_name][read1_ref_name]['rgList'].append(rg)

                    #print(anchor[read2_ref_name][read1_ref_name])

                # CASE 2: read1 maps on a contig and mate read on a chromosome
                elif not 'Chr' in read1_ref_name and 'Chr' in read2_ref_name:

                    read1_posi_list = read.get_reference_positions(full_length=True)
                    read1_start = read1_posi_list[0]
                    read2_start = read.next_reference_start

                    if read1_ref_name not in anchor:
                        anchor[read1_ref_name] = { }

                    if read2_ref_name not in anchor[read1_ref_name]:
                        anchor[read1_ref_name][read2_ref_name] = {
                                                         'posi' :[(read2_start,read1_start)],
                                                         'posiChr': [read2_start],
                                                         'posiCtg': [read1_start],
                                                         'rgList' : [ ],
                                                         'totalHit': 0,
                                                         'startHit': 0,
                                                         'startPosiChr' : [],
                                                         'startPosiCtg': [],
                                                         'endHit': 0,
                                                         'endPosiChr': [],
                                                         'endPosiCtg': [],
                                                         'endPosi' : [],
                                                         'len' :  myRef[read1_ref_name],
                                                         'qual' : []
                                                         }
                    else:
                        anchor[read1_ref_name][read2_ref_name]['posi'].append((read2_start,read1_start))
                        anchor[read1_ref_name][read2_ref_name]['posiChr'].append(read2_start)
                        anchor[read1_ref_name][read2_ref_name]['posiCtg'].append(read1_start)
                        anchor[read1_ref_name][read2_ref_name]['qual'].append(qual)

                    nbNotProp3 += 1
                    anchor[read1_ref_name][read2_ref_name]['totalHit'] += 1
                    if read1_start is not None:
                        if (read1_start < 300):
                            anchor[read1_ref_name][read2_ref_name]['startHit'] += 1
                            anchor[read1_ref_name][read2_ref_name]['startPosiChr'].append(read2_start)
                            anchor[read1_ref_name][read2_ref_name]['startPosiCtg'].append(read1_start)
                            if rg not in  anchor[read1_ref_name][read2_ref_name]['rgList']:
                                anchor[read1_ref_name][read2_ref_name]['rgList'].append(rg)

                        elif (read1_start > myRef[read1_ref_name] - 300):
                            anchor[read1_ref_name][read2_ref_name]['endHit'] += 1
                            anchor[read1_ref_name][read2_ref_name]['endPosiChr'].append(read2_start)
                            anchor[read1_ref_name][read2_ref_name]['endPosiCtg'].append(read1_start)
                            if rg not in  anchor[read1_ref_name][read2_ref_name]['rgList']:
                                anchor[read1_ref_name][read2_ref_name]['rgList'].append(rg)

            # Read is secondary
            if read.is_secondary:
                nbSup += 1

            #read is supplementary
            if read.is_supplementary:
                nbSupTot += 1

            nb += 1

        bamfile.close()
        print(nb, nbNotProp, nbNotProp2, nbNotProp3, nbUnmapped, nbMateUnmapped, nbSup, nbSupTot, nb_noqual)


    ####################################################
    #### PRINT ANCHOR INTO CSV file
    out = open(out_file,"w")
    out.write("CTG_name;CTG_len;CHR_name;ANCHORING_TAG;SAMPLE;SAMPLE_list;TOT_Hit;qual;\
START_Hit;START_CTG_Posi;START_CTG_posi_sd;START_CHR_Posi;START_CHR_posi_sd;START_CHR_Min;\
END_Hit;END_CTG_Posi;END_CTG_posi_sd;END_CHR_Posi;END_CHR_posi_sd;END_CHR_Max;Chr_Ln;\
CHR_all_posi_sd;CHR_all_posi;ALL_POSI_LIST\n")

    for ctg in anchor:
        for chr in anchor[ctg]:

            # To save if the ctg placement on the chr is only in "5'" or "3'" or "5' and 3'"
            tag = ''

            # Get the list of samples from which the reads have been used to place ctgs on the chr
            rgList = anchor[ctg][chr]['rgList']

            # Test that the minimum depth is respectvely in at least on one ctg end
            # According the read depth in the 300first and last pb of the contig, the tag is set
            if anchor[ctg][chr]['startHit'] < depth and anchor[ctg][chr]['endHit'] < depth :
                continue

            elif anchor[ctg][chr]['startHit'] >= depth and anchor[ctg][chr]['endHit'] < depth:
                tag = "5'"

            elif anchor[ctg][chr]['startHit'] < depth and anchor[ctg][chr]['endHit'] >= depth:
                tag = "3'"

            else:
                tag = "5'&3'"

            # Compute the standard deviation on the ctg start positions
            startSdCtg = None

            posiStartCtg = [val for val in anchor[ctg][chr]['startPosiCtg'] if val is not None]
            if posiStartCtg and len(posiStartCtg) > 1:
                startSdCtg = int(nanstd(posiStartCtg))

            # Compute the standard deviation on the chr position (when mate reads maps in the 300st position of the ctg)
            # & chr minimum position
            minPosiChr = None
            startSdChr = None
            startwoOutlierSdChr = None
            startwoOutlierChr = None

            posiStartChr = [val for val in anchor[ctg][chr]['startPosiChr'] if val is not None]
            if posiStartChr :
                # just one chromosome position
                if (len(posiStartChr) == 1):
                    minPosiChr = posiStartChr[0]
                else:
                    # sd on all chromosome positions (with extreme value)
                    startSdChr = int(nanstd(posiStartChr))

                    # To compute sd without extreme values
                    # compute the third quartile and the first quartile.
                    # then extract the positions between these two values
                    low, high = iqr(posiStartChr)
                    startwoOutlierChr = [elem for elem in posiStartChr if elem >= low and elem <= high]

                    # enough values ti compute a standard deviation without extreme ositions?
                    if startwoOutlierChr:
                        startwoOutlierSdChr = int(nanstd(startwoOutlierChr))
                        #print(posiStartChr)
                        #print(startwoOutlierChr)
                        minPosiChr = mean(startwoOutlierChr) # change min
                    else:
                        minPosiChr = mean(posiStartChr)


            # Compute the standard deviation on the ctg end positions
            endSdCtg = None
            posiEndCtg = [val for val in anchor[ctg][chr]['endPosiCtg'] if val is not None]
            if posiEndCtg and len(posiEndCtg) > 1:
                endSdCtg = int(nanstd(posiEndCtg))

            maxPosiChr = None
            endSdChr = None
            endwoOutlierSdChr = None
            endwoOutlierChr = None
            posiEndChr = [val for val in anchor[ctg][chr]['endPosiChr'] if val is not None]
            if posiEndChr :
                if (len(posiEndChr) == 1):
                    maxPosiChr = posiEndChr[0]
                else:
                    endSdChr = int(nanstd(posiEndChr))
                    low, high = iqr(posiEndChr)
                    endwoOutlierChr = [elem for elem in posiEndChr if elem >= low and elem <= high]
                    if endwoOutlierChr:
                        endwoOutlierSdChr = int(nanstd(endwoOutlierChr))
                        #print(posiEndChr)
                        #print(endwoOutlierChr)
                        maxPosiChr = mean(endwoOutlierChr)
                    else:
                        maxPosiChr = mean(posiEndChr)

            posiSdChr = None
            woOutlierSdChr = None
            woOutlierChr = None
            posiChr = [val for val in anchor[ctg][chr]['posiChr'] if val is not None]
            if posiChr:
                ecartChr = int(nanstd(posiChr))
                low, high = iqr(posiChr)
                woOutlierChr = [elem for elem in posiChr if elem >= low and elem <= high]
                if woOutlierChr:
                    woOutlierSdChr = int(nanstd(woOutlierChr))

            chrLen = None
            if tag == "5'&3'" and maxPosiChr is not None and minPosiChr is not None:
                chrLen = maxPosiChr - minPosiChr

            #"CTG name; CTG len; CHR name; ANCHORING TAG; SAMPLE #; SAMPLE list; TOT Hit; \
            #START Hit; START CTG Posi; START CTG posi sd; START CHR Posi ; START CHR posi sd; \
            #END Hit; END CTG Posi; END CTG posi sd; END CHR Posi; END CHR posi sd; \
            #CHR all posi sd; CHR all posi sd; ALL POSI LIST\n")

            line = f"{ctg};{anchor[ctg][chr]['len']};{chr};{tag};{len(rgList)};{ rgList };{anchor[ctg][chr]['totalHit']};{anchor[ctg][chr]['qual']};\
{anchor[ctg][chr]['startHit']};{anchor[ctg][chr]['startPosiCtg']};{startSdCtg };{anchor[ctg][chr]['startPosiChr']};{ startwoOutlierSdChr } ({startSdChr});{ minPosiChr};\
{anchor[ctg][chr]['endHit']};{anchor[ctg][chr]['endPosiCtg']};{endSdCtg};{anchor[ctg][chr]['endPosiChr']};{ endwoOutlierSdChr } ({endSdChr});{maxPosiChr};{chrLen};\
{woOutlierSdChr} ({ ecartChr });{anchor[ctg][chr]['posiChr']};{anchor[ctg][chr]['posi']}\n"
            line = line.replace("None", "-9999")
            line = line.replace(" ", "")

            #print(minPosiChr,maxPosiChr)
            out.write(line)

    out.close()



def main():

    #get arguments
    parser = argparse.ArgumentParser(usage='%(prog)s [options]',
                                     description='Python script to parse bam file ',
                                     epilog="""Author : C. Tranchant-Dubreuil (UMR DIADE / Rice Team - IRD) \n
                                              ,--------------------""")

    parser.add_argument('-b','--bam', help='Bam directory to parse', required=True)
    parser.add_argument('-d','--depth', help='minimum depth', required=True)
    parser.add_argument('-o','--out', help='tabular file', required=True)
    args = parser.parse_args()

    # output file
    out_file = args.out

    # minimum reads depth on a chromosome to place ctgs
    depth = int(args.depth)

    # bam directory with bam and bai files used for ctgs positioning
    bam_dir = args.bam
    file_list = []

    for file in os.listdir(bam_dir):
        full_path = os.path.join(bam_dir, file)
        file_list.append(full_path)

    # Read bam by bam to perform the contig placement according a minimum depth
    read_bam(file_list, depth, out_file)




if __name__ == "__main__":
    main()
