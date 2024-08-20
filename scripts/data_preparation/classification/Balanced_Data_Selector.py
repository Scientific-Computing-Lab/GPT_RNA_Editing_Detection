import sys
import json
import random

def shuffle_jsonl_data(input_file, output_file):
    yes_instances = []
    no_instances = []

    # Read the JSONL file
    with open(input_file, 'r') as file:
        for line in file:
            data = json.loads(line.strip())
            messages = data['messages']
            
            # Extract assistant's response
            assistant_response = messages[-1]['content'].strip()
            
            if assistant_response == "Yes":
                yes_instances.append(data)
            elif assistant_response == "No":
                no_instances.append(data)

    # Randomly select an equal number of "No" instances
    num_yes = len(yes_instances)
    selected_no_instances = random.sample(no_instances, num_yes)

    # Combine and shuffle the selected instances
    combined_instances = yes_instances + selected_no_instances
    random.shuffle(combined_instances)

    # Write shuffled instances to a new JSONL file
    with open(output_file, 'w') as outfile:
        for instance in combined_instances:
            json.dump(instance, outfile)
            outfile.write('\n')

    print(f"Shuffled data with {num_yes} 'Yes' instances and {num_yes} 'No' instances saved to {output_file}.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file.jsonl output_file.jsonl")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    shuffle_jsonl_data(input_file, output_file)
