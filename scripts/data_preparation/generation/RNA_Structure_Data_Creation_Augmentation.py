import sys
from pybedtools import BedTool
import pybedtools
import pandas as pd
import subprocess
import concurrent.futures
import random
import os
import re
from pathlib import Path
from csv import DictReader
import argparse
import csv
import re
from multiprocessing import Pool
from Bio import SeqIO
import json
import alina 
import math
from itertools import combinations, product
import multiprocessing as mp
import argparse
import RNA





class genome_reader:
    """this calss porpose is to enable one time loading of the genome, and use it multiple times to get fasta seqs"""

    def __init__(self, file_p):
        """consturctor - will read fasta file

        Args:
            file_p (str): path to fasta file
        """
        self.data_dict = SeqIO.to_dict(SeqIO.parse(file_p, "fasta"))

    def get_fasta(self, chr, start, end=None):
        """will return fasta seq from start to end if end given,
            otherwise,  will return the the postiaonal nucleotide if only one argument is given for position

        Args:
            chr (str): chr
            start (int): seq wanted start postion if range is given, else - position of wanted nuc
            end (int, optional): end postion if range is given. Defaults to None.

        Returns:
            str: seq
        """
        return (
            str(self.data_dict[chr].seq)[start : end + 1]
            if end
            else str(self.data_dict[chr].seq)[start : start + 1]
        )


def get_reverse_complement(string_seq):
    reverse = string_seq[::-1]
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G','a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'N':'N','n':'n'}
    reverse_complement = ''.join([complement_dict[base] for base in reverse])
    return reverse_complement



def get_seq(chr1,start1,end1,chr2,start2,end2,strand,gr,ds_link):
    side1 = gr.get_fasta(chr1,start1,end1)
    side2 = gr.get_fasta(chr2,start2,end2)
    seq = side1 + ds_link + side2
    seq = seq.upper()
    if strand == "-" :
         seq = get_reverse_complement(seq)
    return seq

def get_adenosines_data(file_path, chr1,start, end):
    adenosines_data = {}
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            position = int(fields[1])
            coverage = int(fields[3])
            chr = fields[0]
            if chr == chr1:
                if start <= position <= end:
                    if coverage > 100:
                        editing_index = float(fields[4])
                        adenosines_data[position] = {'editing_index': round(editing_index, 1), 'coverage': coverage}

    return adenosines_data



def get_data_for_editing_level(tuple_list, index, precent):
    for idx, val in tuple_list:
        if idx == index:
            if val >= precent:
                return "I"
            else:
                return "A"
    return "A"  


def is_within_edit_range(base_pair, edit_site, range_size):
    """
    Function to check if a base pair is within a specified range around an edit site.
    """
    edit_index, edit_level = edit_site
    base_index, base_comp_index = base_pair
    return abs(edit_index - base_index) <= range_size or abs(edit_index - base_comp_index) <= range_size

def filter_base_pairs(base_pairs, edit_sites, range_size):
    """
    Function to filter base pairs based on edit site proximity.
    """
    filtered_pairs = []
    for base_pair in base_pairs:
        within_range = False
        for edit_site in edit_sites:
            if is_within_edit_range(base_pair, edit_site, range_size):
                within_range = True
                break
        if not within_range:
            filtered_pairs.append(base_pair)
    return filtered_pairs


def calculate_identity_percentage(string1, string2):
    if len(string1) != len(string2):
        raise ValueError("Strings must have the same length")
    
    # Count matching characters
    match_count = sum(c1 == c2 for c1, c2 in zip(string1, string2))
    
    # Calculate percentage of identity
    percentage_identity = (match_count / len(string1)) * 100
    
    return percentage_identity


         
def augment_data_func(seq_info, paired_bases, editing_sites, editing_level,distance,percentage_change, method,percent_structural_identity):
    vienna = seq_info[2]
    new_seq = list(seq_info[1])

    # Filter editing sites that meet the threshold
    eligible_editing_sites = [(site, level) for site, level in editing_sites if level >= editing_level]

    # Calculate maximum number of allowed changes
    allowed_pair = filter_base_pairs(paired_bases, eligible_editing_sites, distance)
    num_allowed_pair = len(allowed_pair)
    num_changes = math.ceil((percentage_change * num_allowed_pair) / 100 )

    change_position_idx_list = []

    for i in range(num_changes):
        # Select position for change
        if method == 'random':
            change_position_idx = random.randint(0, num_allowed_pair-1)
            while change_position_idx  in change_position_idx_list:
                change_position_idx = random.randint(0, num_allowed_pair-1)
            change_position_idx_list.append(change_position_idx)
            change_position = allowed_pair[change_position_idx]
        else:
            allowed_pair_sort = sorted(allowed_pair, key=lambda x: x[0])
            change_position = allowed_pair_sort[i]
        pos = change_position[0]
        comp_pos = change_position[1]

        if seq_info[1][pos] == 'G':
            new_seq[pos] = 'C'
            new_seq[comp_pos] = 'G'
        elif seq_info[1][pos] == 'C':
            new_seq[pos] = 'G'
            new_seq[comp_pos] = 'C'
        elif seq_info[1][pos] == 'T':
            new_seq[pos] = 'A'
            new_seq[comp_pos] = 'T'
        else:
            new_seq[pos] = 'T'
            new_seq[comp_pos] = 'A'

    new_seq_s = ''.join(new_seq)

    fc = RNA.fold_compound(new_seq_s)
    # compute MFE and MFE structure
    (mfe_struct, mfe) = fc.mfe()
    dbn_data = (mfe,new_seq_s,mfe_struct)
    percentage_identity = calculate_identity_percentage(dbn_data[2], vienna)
    return  dbn_data, percentage_identity

def process_row( row,gr,A_pos_A2G_file,A_pos_T2C_file,editing_level,distance,percent_structural_identity,min_percent_change,max_percent_change,num_augment_seq,method):
    chr1, start1, end1, chr2, start2, end2, strand = row

    ds_link = "NNNNNNNNNNNNNNN"
    seq = get_seq(chr1, start1, end1, chr2, start2, end2, strand, gr,ds_link)
    
    fc = RNA.fold_compound(seq)
    # compute MFE and MFE structure
    (mfe_struct, mfe) = fc.mfe()
    dbn_data = (mfe,seq,mfe_struct)


    if strand == '+':
        adenosines_data = get_adenosines_data(A_pos_A2G_file, chr1,start1, end2)
        alu1_data_for_output = [(position - start1 + 1,adenosines_data[position]['editing_index']) for position in adenosines_data if start1 <= position <= end1]
        alu2_data_for_output = [((end1-start1) + position - start2 + len(ds_link) + 2,adenosines_data[position]['editing_index']) for position in adenosines_data if start2 <= position <= end2]

    else:
        adenosines_data = get_adenosines_data(A_pos_T2C_file, chr1,start1, end2)
        alu2_data_for_output = [(end2- position + 1,adenosines_data[position]['editing_index']) for position in adenosines_data if start2 <= position <= end2]
        alu1_data_for_output = [( (end2-start2) + end1 - position + len(ds_link) + 2, adenosines_data[position]['editing_index']) for position in adenosines_data if start1 <= position <= end1]

    
    all_A_data = alu1_data_for_output + alu2_data_for_output

    struct=dbn_data[2]
    pairs = alina.utils.get_na_pairs(struct)
    predicate0 = "'"
    for i, base in enumerate(dbn_data[1]):
        if base != "A":
            predicate0 += base
        else:
            predicate0 += get_data_for_editing_level(all_A_data,i+1,editing_level)
    predicate0 += "'"
    augment_data = [[chr1, start1, end1, chr2, start2, end2,strand, dbn_data[1:], predicate0,0,100]]

    for i in range(num_augment_seq):
        percentage_change = random.randint(min_percent_change, max_percent_change)
        augment_dbn_data, percentage_identity = augment_data_func(dbn_data, pairs, all_A_data, editing_level, distance, percentage_change, method, percent_structural_identity)
        while percentage_identity < percent_structural_identity:
            augment_dbn_data, percentage_identity = augment_data_func(dbn_data, pairs, all_A_data, editing_level, distance, percentage_change, method, percent_structural_identity, fr)
        predicate_editing_site = "'"
        for i, base in enumerate(augment_dbn_data[1]):
            if base != "A":
                predicate_editing_site += base
            else:
                predicate_editing_site += get_data_for_editing_level(all_A_data,i+1,editing_level)
        predicate_editing_site += "'"

        augment_data.append([chr1, start1, end1, chr2, start2, end2,strand, augment_dbn_data[1:], predicate_editing_site,percentage_change, percentage_identity])
    return augment_data

if __name__ == "__main__":

    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(formatter_class=MyFormatter,
                                     description='Creating a dataset of secondary structures of editing sites  and augmentation.')
    parser.add_argument('--pair_region', '--pair_region', dest="pair_region", help='A bed file that contains both alu (chr1,start1,end1,chr2,start2,end2,strand)', type=str,required=True)
    parser.add_argument('--output_dir', '--output_dir', dest="output_dir", help='dir for the output file (csv)', type=str,required=True)
    parser.add_argument('--editing_site_plus', '--editing_site_plus', dest="editing_site_plus", help='A file that contains data on the levels of editing in adenosines in the positive strand', type=str,required=True)
    parser.add_argument('--editing_site_minus', '--editing_site_minus', dest="editing_site_minus", help='A file that contains data on the levels of editing in adenosines in the negative strand', type=str,required=True)
    parser.add_argument('--genome', '--genome', dest="genome", help='fasta file of the genome',type=str, required=True)
    parser.add_argument('--processes', '--processes', dest="processes", help='number of processes', required=True,type=int)
    parser.add_argument('--editing_level', '--editing_level', dest="editing_level", help='The minimum editing level for a site to be considered an editing site', required=True,type=float)
    parser.add_argument('--distance', '--distance', dest="distance", help='The distance from an editing site where there will be no sequential change in augmentation', required=True,type=int)
    parser.add_argument('--num_augment_seq', '--num_augment_seq', dest="num_augment_seq", help='The number of new sequences that will be created from each original sequence', required=True,type=int)
    parser.add_argument('--min_percent_change', '--min_percent_change', dest="min_percent_change", help='The minimum percentage of the structure sequence in which there will be a change in augmentation', required=True,type=int)
    parser.add_argument('--max_percent_change', '--max_percent_change', dest="max_percent_change", help='The maximum percentage of the structure sequence in which there will be a change in augmentation', required=True,type=int)
    parser.add_argument('--structural_identity', '--structural_identity', dest="percent_structural_identity", help='The minimum identity percentage in Vienna for an augmented structure', required=True,type=int)
    parser.add_argument('--method', '--method', dest="method", choices = ['random','linear'], help='The method of selecting the change sites', required=True, metavar='method')

    args = parser.parse_args()

    alu_region = pd.read_csv(args.pair_region, sep='\t', header=None)

    alu_region.columns = ['chr1', 'start1', 'end1','chr2','start2','end2','strand']

    gr = genome_reader(args.genome)
    list_for_row = alu_region.values.tolist()


    with Pool(processes=args.processes) as pool:
        results = pool.starmap(process_row,[(list_for_row[i],gr,args.editing_site_plus,args.editing_site_minus,args.editing_level,args.distance,args.percent_structural_identity,args.min_percent_change,args.max_percent_change,args.num_augment_seq,args.method) for i in range(len(list_for_row))])

    pool.close()
    pool.join()
    
    all_augment_data = [item for sublist in results for item in sublist]
    df_results_augment_seq = pd.DataFrame(all_augment_data)
    df_final_result_path = os.path.join(args.output_dir + "all_final_result.csv")
    df_results_augment_seq.to_csv(df_final_result_path, index=False, header=False)
    df_for_train = df_results_augment_seq.iloc[:,7:9]
    data_for_prepare = os.path.join(args.output_dir + "data_for_prepare.csv")
    df_for_train.to_csv(data_for_prepare, index=False, header=['vienna_format', 'feature_annotation'])

        


