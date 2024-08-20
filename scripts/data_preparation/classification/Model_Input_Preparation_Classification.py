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
        "content": "Predict if the central adenosine (A) in the given RNA sequence context within an Alu element will be edited to inosine (I) by ADAR enzymes."
    }

    data = []

    with open(input_filename, mode='r', newline='') as csv_file:
        csv_reader = csv.reader(csv_file)
        header = next(csv_reader)  # Skip the header row if there is one

        for row in csv_reader:
            structure,L,R,y_n = row
            contant = f"L:" + L + ", A:A, R:" + R + ", Alu Vienna Structure:" + structure 
            data.append({
                "messages": [
                    system_message,
                    {"role": "user", "content": contant},
                    {"role": "assistant", "content": y_n}
                ]
            })

    # Shuffle the data
    random.shuffle(data)

    # Split the data into 80% training and 20% evaluation
    split_index = int(len(data) * 0.8)
    training_data = data[:split_index]
    evaluating_data = data[split_index:]
    
    # Write the training data to the training JSONL file
    with open(training_filename, mode='w') as train_file:
        for entry in training_data:
            train_file.write(json.dumps(entry) + '\n')

    # Write the evaluating data to the evaluating JSONL file
    with open(evaluating_filename, mode='w') as eval_file:
        for entry in evaluating_data:
            eval_file.write(json.dumps(entry) + '\n')

    print(f"Training data saved to {training_filename}.")
    print(f"Evaluation data saved to {evaluating_filename}.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python csv_to_jsonl.py <input_csv_file>")
        sys.exit(1)

    input_filename = sys.argv[1]
    csv_to_jsonl(input_filename)

