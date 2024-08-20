import os
import sys
import json
import re
from openai import AzureOpenAI

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 4:
    print("Usage: python script.py <input_file_path> <output_file_path> <temperature>")
    sys.exit(1)

input_file_path = sys.argv[1]
output_file_path = sys.argv[2]
temperature_test = float(sys.argv[3])

# Initialize the Azure OpenAI client
client = AzureOpenAI(
    azure_endpoint="https://ncusgaloren.openai.azure.com/",
    api_key=os.getenv("AZURE_OPENAI_KEY"),
    api_version="2024-02-15-preview"
)

# Initialize meta-statistics counts
total_fp = total_fn = total_tp = total_tn = 0

# Create a log file
log_file = os.path.splitext(output_file_path)[0] + "_log.txt"

def log(message):
    with open(log_file, 'a') as logf:
        logf.write(message + '\n')

# Function to compare two Alu strings and calculate FP, FN, TP, TN
def compare_annotation_strings(true_string, predicted_string):
    fp = fn = tp = tn = 0
    
    for true_char, predicted_char in zip(true_string, predicted_string):
        if true_char == 'I' and predicted_char == 'I':
            tp += 1
        elif true_char != 'I' and predicted_char == 'I':
            fp += 1
        elif true_char == 'I' and predicted_char != 'I':
            fn += 1
        else:
            tn += 1
    
    return fp, fn, tp, tn

# Function to process each line of the JSONL file, generate predictions, and compare results
def extract_and_process_jsonl(input_path, output_path):
    global total_fp, total_fn, total_tp, total_tn
    
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            data = json.loads(line)
            messages = data.get('messages', [])
            user_content = ""
            assistant_content = ""

            for message in messages:
                if message['role'] == 'user':
                    user_content = message['content']
                elif message['role'] == 'assistant':
                    assistant_content = message['content']

            if user_content and assistant_content:
                log("Processing user content: " + user_content)
                
                # Constructing message_text as a list of dictionaries
                message_text = [
                    {
                        "role": "system",
                        "content": "Your mission is to assist in the analysis of RNA Alu sequences and the identification of ADAR (A-to-I) editing locations."
                    },
                    {
                        "role": "user",
                        "content": user_content
                    }
                ]
                
                log("Sending request to OpenAI API...")
                
                response = client.chat.completions.create(
                    model="15editingLevel_cov100",
                    messages=message_text,
                    temperature=temperature_test,
                    max_tokens=4096,
                    logprobs=True,
                    top_logprobs=5,
                    top_p=0.95,
                    frequency_penalty=0,
                    presence_penalty=0,
                    stop=None
                )

                log("Received response from OpenAI API.")
                print(response)
                generated_response = response.choices[0].message.content
                
                log("Generated response: " + generated_response)
                log("Comparing true and predicted strings...")
                
                fp, fn, tp, tn = compare_annotation_strings(assistant_content, generated_response)
                
                log("*******************************")
                log(f"FP, FN, TP, TN: {fp}, {fn}, {tp}, {tn}")
                log("*******************************")
                
                # Update meta-statistics counts
                total_fp += fp
                total_fn += fn
                total_tp += tp
                total_tn += tn

                result = {
                    "user_content": user_content,
                    "assistant_content": assistant_content,
                    "generated_response": generated_response,
                    "fp": fp,
                    "fn": fn,
                    "tp": tp,
                    "tn": tn
                }

                outfile.write(json.dumps(result) + '\n')
                log("Result saved to output file.")
                log("-" * 50)

# Run the function to process the JSONL file
log("Starting processing...")
extract_and_process_jsonl(input_file_path, output_file_path)

# Calculate meta-statistics
total_samples = total_fp + total_fn + total_tp + total_tn
meta_statistics = {
    "total_samples": total_samples,
    "total_fp": total_fp,
    "total_fn": total_fn,
    "total_tp": total_tp,
    "total_tn": total_tn,
    "meta_fp": total_fp / total_samples,
    "meta_fn": total_fn / total_samples,
    "meta_tp": total_tp / total_samples,
    "meta_tn": total_tn / total_samples
}

log("\nMeta-Statistics:")
log(json.dumps(meta_statistics, indent=4))

# Save meta-statistics in a separate file
meta_stats_file = os.path.splitext(output_file_path)[0] + "_meta_stats.json"
with open(meta_stats_file, 'w') as metafile:
    json.dump(meta_statistics, metafile)

log(f"Meta-statistics saved in {meta_stats_file}")
