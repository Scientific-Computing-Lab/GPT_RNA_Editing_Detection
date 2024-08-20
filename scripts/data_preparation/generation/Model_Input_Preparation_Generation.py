import csv
import json
import random
import os
import sys
from datetime import datetime

def csv_to_jsonl(input_filename):
    base_name = os.path.splitext(input_filename)[0]
    
    # Generate a unique identifier using the current date and time
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Create filenames with the unique timestamp
    output_filename = f"{base_name}_{timestamp}.jsonl"
    training_filename = f"train_{base_name}_{timestamp}.jsonl"
    evaluating_filename = f"valid_{base_name}_{timestamp}.jsonl"

    system_message = {
        "role": "system",
        "content": "Your mission is to assist in the analysis of RNA Alu sequences and the identification of ADAR (A-to-I) editing locations."
    }

    data = []

    with open(input_filename, mode='r', newline='') as csv_file:
        csv_reader = csv.reader(csv_file)
        header = next(csv_reader)  # Skip the header row if there is one

        for row in csv_reader:
            if len(row) != 2:
                continue  # Skip rows that do not have exactly two columns

            nucleotide_sequence, feature_annotation = row
            data.append({
                "messages": [
                    system_message,
                    {"role": "user", "content": nucleotide_sequence},
                    {"role": "assistant", "content": feature_annotation}
                ]
            })
 

    # Split the data into 80% training and 20% evaluation
    split_index = int(len(data) * 0.8)
    training_data = data[:split_index]
    evaluating_data = data[split_index:]
    
    # Calculate the number of samples for 10% of the training and evaluating data
    num_train_samples = int(len(training_data) * 0.1)
    num_eval_samples = int(len(evaluating_data) * 0.1)

    # Randomly select 10% of the samples from training_data and evaluating_data
    random_train_samples = random.sample(training_data, num_train_samples)
    random_eval_samples = random.sample(evaluating_data, num_eval_samples)

    # Write the training data to the training JSONL file
    with open(training_filename, mode='w') as train_file:
        for entry in random_train_samples:
            train_file.write(json.dumps(entry) + '\n')

    # Write the evaluating data to the evaluating JSONL file
    with open(evaluating_filename, mode='w') as eval_file:
        for entry in random_eval_samples:
            eval_file.write(json.dumps(entry) + '\n')

    print(f"Training data saved to {training_filename}.")
    print(f"Evaluation data saved to {evaluating_filename}.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python prepare.py <input_csv_file>")
        sys.exit(1)

    input_filename = sys.argv[1]
    csv_to_jsonl(input_filename)

