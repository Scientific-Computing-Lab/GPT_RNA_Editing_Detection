import json
#pip install ViennaRNA
import RNA
import sys




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

def find_min_max(numbers):
    if not numbers:
        return None, None  # Return None for both min and max if the list is empty
    return min(numbers), max(numbers)


def identify_loops(sequence, structure):
    paired_bases = parse_secondary_structure(sequence, structure)
    loops = []
    
    i = 1  # RNA.ptable uses 1-based indexing
    while i <= len(sequence):
        if paired_bases[i] == 0:  # Unpaired base (part of a loop)
            loop = []
            while i <= len(sequence) and paired_bases[i] == 0:
                loop.append(i - 1)  # Convert to 0-based indexing
                i += 1
            loops.append(loop)
        else:
            i += 1
    
    return loops

def find_loop_pairing(sequence, structure):
    loops = identify_loops(sequence, structure)
    loop_pairs = []
    
    for loop in loops:
        loop_info = {'loop_bases': [],'loop_index' : [], 'loop_pairs': [], 'pair_index' : []}
        
        # Collect loop bases
        for i in loop:
            loop_info['loop_bases'].append(sequence[i])
            loop_info['loop_index'].append(i)
        min_number, max_number = find_min_max(loop)
        min_number = min_number -1
        max_number = max_number + 1
        paired_bases_index = parse_dot_bracket(structure)
        min_index = paired_bases_index[min_number]
        max_index = paired_bases_index[max_number]
        paired_index = []
        paired_seq = []
        for i in range(max_index+1,min_index):
            paired_index.append(i)
            paired_seq.append(sequence[i])
            loop_info['pair_index'].append(i)
            loop_info['loop_pairs'].append(sequence[i])
        if len(paired_seq) > 1: 
            loop_info['loop_pairs'].reverse() 
            loop_info['pair_index'].reverse()
            paired_seq.reverse()
            paired_index.reverse()       
        loop_pairs.append(loop_info)
    
    return loop_pairs


# Function to parse the secondary structure and determine the pairing status of each base
def parse_secondary_structure(sequence, structure):
    paired_bases = RNA.ptable(structure)
    return paired_bases

# Function to check the editing site conditions
def check_editing_site(structure, L, R):
    # Combine L, A, R to get the full sequence and find the position of A
    full_sequence = L + "A" + R
    central_A_index = len(L)
    # Parse the secondary structure
    loop_pairs_info = find_loop_pairing(full_sequence, structure)
    # Check if the central A is in a loop
    befor_A_seq = full_sequence[central_A_index-1]
    pair_base = ''
    for d in loop_pairs_info:
        loop_index_list = d.get('loop_index', [])  # Get the 'loop_index' list from the dictionary
        pair_list = d.get('loop_pairs', [])
        for idx, num in enumerate(loop_index_list):
            if num == central_A_index-1:
                if len(pair_list) == 0:
                    pair_base = ''
                    break
                if len(pair_list) != len(loop_index_list):
                    pair_base = ''
                else:
                    pair_base = pair_list[idx]
                break
    
    if pair_base == 'C' and befor_A_seq != 'G':
        return True

    return False


if __name__ == "__main__":

    TP = 0  # True Positives
    FP = 0  # False Positives
    TN = 0  # True Negatives
    FN = 0  # False Negatives

    input_filename = sys.argv[1]
    # Load the JSON Lines data
    with open(input_filename, 'r') as file:
        lines = file.readlines()

    # Process each JSON object
    for i, line in enumerate(lines):
        data = json.loads(line.strip())

        # Extract the relevant information from the JSON data
        for message in data["messages"]:
            assistant = data["messages"][2] 
            Editing_site = assistant["content"]
            if message["role"] == "user":
                content = message["content"]
                # Split the content to get sequences and Vienna structure
                parts = content.split(", ")
                L = parts[0].split(":")[-1].strip()
                R = parts[2].split(":")[-1].strip()
                vienna_structure = parts[-1].split("Alu Vienna Structure:")[-1].strip()

                # Check if the sequence meets the editing site conditions
                is_potential_editing_site = check_editing_site(vienna_structure, L, R)


                if Editing_site == "Yes" and is_potential_editing_site:
                    TP +=1
                if Editing_site == "Yes" and not is_potential_editing_site:
                    FN+=1
                if Editing_site == "No" and is_potential_editing_site:
                    FP+=1
                if Editing_site == "No" and not is_potential_editing_site:
                    TN +=1


    accuracy = ((TP + TN) / (TP+TN+FP+FN)) * 100
    specifity = (TN / (TN+FP)) * 100
    precision = (TP / (TP+FP)) * 100
    recall = (TP / (TP+FN)) * 100
    f1score = ((2 * (recall*precision))/ (recall + precision))

    print("Total TP:" + str(TP))
    print("Total FN:" + str(FN))
    print("Total FP:" + str(FP))
    print("Total TN:" + str(TN))
    print("Accuracy:" + str(round(accuracy,2)) + "%")
    print("Precision:" + str(round(precision,2)) + "%")
    print("Recall:" + str(round(recall,2)) + "%")
    print("Specificity:" + str(round(specifity,2)) + "%")
    print("F1 Score:" + str(round(f1score,2)) + "%")


