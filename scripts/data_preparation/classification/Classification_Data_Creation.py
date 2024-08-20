import pandas as pd
import subprocess
import argparse
import csv
from multiprocessing import Pool
import argparse
import pandas as pd
import subprocess
import argparse
import csv
from multiprocessing import Pool
import multiprocessing as mp
import argparse
import RNA
import os
from Bio import SeqIO



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



def get_seq(chr1,start1,end1,chr2,start2,end2,strand,gr):
    side1 = gr.get_fasta(chr1,start1,end1)
    side2 = gr.get_fasta(chr2,start2,end2)
    seq = side1 + "NNNNNNNNNN" + side2
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



def count_letter(sequence, letter):
    count = 0
    for char in sequence:
        if char == letter:
            count += 1
    return count

def get_data_for_editing_level(tuple_list, index,editing_level):
    for idx, val in tuple_list:
        if idx == index:
            if val >= float(editing_level):
                return "Yes"
            else:
                return "No"
    return "No"

def parse_dot_bracket(structure):
    """
    Parses dot-bracket notation to identify base pairs.
    
    Args:
    - structure (str): The RNA secondary structure in dot-bracket notation.
    
    Returns:
    - dict: A dictionary with base pair indices.
    """
    stack = []
    pairs = {}
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs[i] = j
                pairs[j] = i
    return pairs

def find_adenosines(sequence):
    """
    Find all the positions of adenosines (A) in the sequence.
    
    Args:
    - sequence (str): The RNA sequence.
    
    Returns:
    - List[int]: A list of positions where 'A' is found.
    """
    return [i for i, nucleotide in enumerate(sequence) if nucleotide == 'A']

def get_subranges(index, length,ds,first_N,last_N, window=40):
    """
    Gets the subranges around a given index.
    
    Args:
    - index (int): The position of the nucleotide (adenosine).
    - length (int): Length of the sequence.
    - window (int): Number of nucleotides to include on either side of the index.
    
    Returns:
    - List[int]: The list of indices around the given index.
    """
    if ds == 1:
        start = max(0, index - window)
        end = min(first_N, index + window + 1)
    else:
        start = max(last_N, index - window)
        end = min(length, index + window + 1)
    return list(range(start, end))

def extract_substructure(sequence, structure, indices):
    """
    Extracts the subsequence and structure for the given indices.
    
    Args:
    - sequence (str): The RNA sequence.
    - structure (str): The RNA secondary structure in dot-bracket notation.
    - indices (List[int]): The list of indices to extract.
    
    Returns:
    - (str, str): The subsequence and its corresponding secondary structure.
    """
    sub_sequence = ''.join(sequence[i] for i in indices)
    sub_structure = ''.join(structure[i] for i in indices)
    
    return sub_sequence, sub_structure

def find_first_last_N(s):
    first_index = s.find('N')
    last_index = s.rfind('N')
    return first_index, last_index


def check_vienna_dot_bracket(sequence, structure):
    # Check if the lengths of the sequence and structure match
    if len(sequence) != len(structure):
        return False

    # Check if the dot-bracket notation is balanced
    stack = []
    for char in structure:
        if char == '(':
            stack.append('(')
        elif char == ')':
            if not stack or stack.pop() != '(':
                return False

    if stack:
        return False

    return True

def check_no_overlap(list1, list2):
    # Find the minimum and maximum of each list
    min1, max1 = list1[0], list1[-1]
    min2, max2 = list2[0], list2[-1]

    # Check if there is no overlap
    if max1 < min2 or max2 < min1:
        return True
    else:
        return False
    
def generate_full_context(sequence, structure,all_A_data,chr1,start1,end1,chr2,start2,end2,strand,editing_level):
    """
    Generates the full_context Vienna format for adenosines in the sequence.
    
    """
    adenosine_positions = find_adenosines(sequence)
    full_pairs = parse_dot_bracket(structure)
    substructures = []
    first_N, last_N = find_first_last_N(sequence)
    ds = 0
    for position in adenosine_positions:
        if position < first_N:
            ds = 1
        else:
            ds = 2
        L = sequence[0:position]
        R = sequence[position+1:]
        editing_level_y_n = get_data_for_editing_level(all_A_data,position+1,editing_level)
        substructures.append((chr1,start1,end1,chr2,start2,end2,strand,sequence,0,0, structure,L,R,editing_level_y_n,position))

    
    return substructures


def generate_abbreviated_vienna(sequence, structure,all_A_data,chr1,start1,end1,chr2,start2,end2,strand,editing_level,ds_link):
    """
    Generates the abbreviated Vienna format for adenosines in the sequence.

    """
    adenosine_positions = find_adenosines(sequence)
    full_pairs = parse_dot_bracket(structure)
    substructures = []
    first_N, last_N = find_first_last_N(sequence)
    ds = 0
    for position in adenosine_positions:
        if position == 0:
            continue
        if position < first_N:
            ds = 1
        else:
            ds = 2
        short = False
        subrange1 = get_subranges(position, len(sequence),ds,first_N,last_N)
        if subrange1[0] not in full_pairs or subrange1[len(subrange1)-1] not in full_pairs:
            continue
       
        opposite_indices = []
        for i in subrange1:
            if i in full_pairs:
                opposite_indices.append(full_pairs[i])
        
        if opposite_indices:
            start_opposite = min(opposite_indices)
            end_opposite = max(opposite_indices) + 1
            subrange2 = list(range(start_opposite, end_opposite))
        else:
            subrange2 = []
        if subrange1[0] > subrange2[0]:
            s = subrange1
            subrange1 = subrange2
            subrange2 = s
        sub_seq1, sub_struct1 = extract_substructure(sequence, structure, subrange1)
        sub_seq2, sub_struct2 = extract_substructure(sequence, structure, subrange2)
        
        # Connect the two substructures with 'NNNNNNNN'
        combined_sequence = sub_seq1 + ds_link + sub_seq2
        combined_structure = sub_struct1 + '..........' + sub_struct2
        is_correct = check_vienna_dot_bracket(combined_sequence, combined_structure)
        no_overlap = check_no_overlap(subrange1, subrange2)

        if not is_correct or not no_overlap:
            continue
        if ds == 1:
            up_pos = position - subrange1[0]
            L = combined_sequence[0:up_pos]
            R = combined_sequence[up_pos+1:]
        else:
            up_pos = position - subrange2[0]
            L = combined_sequence[0:len(sub_seq1) + len(ds_link) + up_pos]
            R = combined_sequence[len(sub_seq1) + len(ds_link) + up_pos+1:]
        editing_level_y_n = get_data_for_editing_level(all_A_data,position+1,editing_level)
        substructures.append((chr1,start1,end1,chr2,start2,end2,strand,sequence, structure,combined_sequence,combined_structure,L,R,editing_level_y_n,position,sub_seq1,sub_struct1,sub_seq2,sub_struct2))
    
    return substructures


if __name__ == '__main__':

    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass


    parser = argparse.ArgumentParser(formatter_class=MyFormatter,
                                     description='Creating a dataset of secondary structures of editing sites  and augmentation.')
    parser.add_argument('--pair_region', '--pair_region', dest="pair_region", help='A bed file that contains both alu (chr1,start1,end1,chr2,start2,end2,strand)', type=str,required=True)
    parser.add_argument('--output_dir', '--output_dir', dest="output_dir", help='dir for the output file (csv)', type=str,required=True)
    parser.add_argument('--editing_site_plus', '--editing_site_plus', dest="editing_site_plus", help='A file that contains data on the levels of editing in adenosines in the positive strand', type=str,required=True)
    parser.add_argument('--editing_site_minus', '--editing_site_minus', dest="editing_site_minus", help='A file that contains data on the levels of editing in adenosines in the negative strand', type=str,required=True)
    parser.add_argument('--genome', '--genome', dest="genome", help='fasta file of the genome',type=str, required=True)
    parser.add_argument('--editing_level', '--editing_level', dest="editing_level", help='The minimum editing level for a site to be considered an editing site', required=True,type=float)
    parser.add_argument('--context', '--context', dest="context", choices = ['full','partial'], help='partial/full context ', required=True, metavar='context')

    args = parser.parse_args()

    alu_region = pd.read_csv(args.pair_region, sep='\t', header=None)

    alu_region.columns = ['chr1', 'start1', 'end1','chr2','start2','end2','strand']

    gr = genome_reader(args.genome)

    all_sub_alu = []
    ds_link = "NNNNNNNNNN"

    for index, row in alu_region.iterrows():
        chr1, start1, end1, chr2, start2, end2, strand = row

        seq = get_seq(chr1,start1,end1,chr2,start2,end2,strand,gr)
        fc = RNA.fold_compound(seq)
        # compute MFE and MFE structure
        (mfe_struct, mfe) = fc.mfe()
        dbn_data = (mfe,seq,mfe_struct)


        if strand == '+':
            adenosines_data = get_adenosines_data(args.editing_site_plus, chr1,start1, end2)
            alu1_data = [(position, adenosines_data[position]['editing_index'], position - start1 + 1, adenosines_data[position]['coverage']) for position in adenosines_data if start1 <= position <= end1]
            alu1_data_for_output = [(position - start1 + 1,adenosines_data[position]['editing_index']) for position in adenosines_data if start1 <= position <= end1]

            alu2_data = [(position, adenosines_data[position]['editing_index'], (end1-start1) + position - start2 + len(ds_link) + 2, adenosines_data[position]['coverage']) for position in adenosines_data if start2 <= position <= end2]
            alu2_data_for_output = [((end1-start1) + position - start2 + len(ds_link) + 2,adenosines_data[position]['editing_index']) for position in adenosines_data if start2 <= position <= end2]

        else:
            adenosines_data = get_adenosines_data(args.editing_site_minus, chr1,start1, end2)
            alu2_data = [(position, adenosines_data[position]['editing_index'], end2- position + 1, adenosines_data[position]['coverage']) for position in adenosines_data if start2 <= position <= end2]
            alu2_data_for_output = [(end2- position + 1,adenosines_data[position]['editing_index']) for position in adenosines_data if start2 <= position <= end2]

            alu1_data = [(position, adenosines_data[position]['editing_index'], (end2-start2) + end1 - position + len(ds_link) + 2, adenosines_data[position]['coverage']) for position in adenosines_data if start1 <= position <= end1]
            alu1_data_for_output = [( (end2-start2) + end1 - position + len(ds_link) + 2, adenosines_data[position]['editing_index']) for position in adenosines_data if start1 <= position <= end1]

    
        all_A_data = alu1_data_for_output + alu2_data_for_output
        sequence = dbn_data[1]
        structure = dbn_data[2]

    # Generate abbreviated Vienna formats
        if args.context == "partial":
            abbreviated_viennas = generate_abbreviated_vienna(sequence, structure,all_A_data,chr1,start1,end1,chr2,start2,end2,strand,args.editing_level,ds_link)
            all_sub_alu = all_sub_alu + abbreviated_viennas
        if args.context == "full":
            abbreviated_viennas = generate_full_context(sequence, structure,all_A_data,chr1,start1,end1,chr2,start2,end2,strand,args.editing_level)
            all_sub_alu = all_sub_alu + abbreviated_viennas           

    df_results = pd.DataFrame(all_sub_alu)
    df_final_result_path = os.path.join(args.output_dir + "all_data_classification_" + args.context + ".csv")
    df_results.to_csv(df_final_result_path, index=False, header=False)
    df_for_train = df_results.iloc[:,10:14]
    data_for_prepare = os.path.join(args.output_dir + "data_for_prepare_classification_" + args.context + ".csv")
    df_for_train.to_csv(data_for_prepare, index=False, header=['structure','L','R','editing_level'])





