import os
import sys
import json
import re
import math
from openai import AzureOpenAI
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score, confusion_matrix
import numpy as np

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

# Initialize meta-statistics counts and probability sums
total_fp = total_fn = total_tp = total_tn = 0
sum_prob_fp = sum_prob_fn = sum_prob_tp = sum_prob_tn = 0
true_labels = []
predicted_labels = []
probabilities = []

# Create a log file
log_file = os.path.splitext(output_file_path)[0] + "_log.txt"

def log(message):
    with open(log_file, 'a') as logf:
        logf.write(message + '\n')

# Function to compare two annotation strings and calculate FP, FN, TP, TN
def compare_annotation_strings(true_label, predicted_label, prob_yes):
    if true_label == "Yes" and predicted_label == "Yes":
        return 0, 0, 1, 0, 0, 0, prob_yes, 0
    elif true_label == "No" and predicted_label == "Yes":
        return 1, 0, 0, 0, prob_yes, 0, 0, 0
    elif true_label == "Yes" and predicted_label == "No":
        return 0, 1, 0, 0, 0, prob_yes, 0, 0
    else:
        return 0, 0, 0, 1, 0, 0, 0, prob_yes

# Function to process each line of the JSONL file, generate predictions, and compare results
def extract_and_process_jsonl(input_path, output_path):
    global total_fp, total_fn, total_tp, total_tn
    global sum_prob_fp, sum_prob_fn, sum_prob_tp, sum_prob_tn
    global true_labels, predicted_labels, probabilities
    
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
                        "content": "Predict if the central adenosine (A) in the given RNA sequence context within an Alu element will be edited to inosine (I) by ADAR enzymes."
                    },
                    {
                        "role": "user",
                        "content": user_content
                    }
                ]
                
                log("Sending request to OpenAI API...")
                
                response = client.chat.completions.create(
                    model="15_yesno_def_EQ",
                    messages=message_text,
                    temperature=temperature_test,
                    max_tokens=1,
                    logprobs=True,
                    top_logprobs=2,
                    top_p=1,
                    frequency_penalty=0,
                    presence_penalty=0,
                    stop=None
                )

                log("Received response from OpenAI API.")
                generated_response = response.choices[0].message.content
                
                top_logprobs = response.choices[0].logprobs.content[0].top_logprobs
                logprob_yes = next((prob.logprob for prob in top_logprobs if prob.token == 'Yes'), None)
                logprob_no = next((prob.logprob for prob in top_logprobs if prob.token == 'No'), None)

                if logprob_yes is not None:
                    prob_yes = math.exp(logprob_yes)
                else:
                    prob_yes = 0.0
                
                if logprob_no is not None:
                    prob_no = math.exp(logprob_no)
                else:
                    prob_no = 0.0

                log(f"Generated response: {generated_response}")
                #log(f"Top Logprobs - Yes: {logprob_yes}, No: {logprob_no}")
                
                log("Comparing true and predicted strings...")
                
                fp, fn, tp, tn, prob_fp, prob_fn, prob_tp, prob_tn = compare_annotation_strings(assistant_content, generated_response, prob_yes)
                
                log("*******************************")
                log(f"FP, FN, TP, TN: {fp}, {fn}, {tp}, {tn}")
                log(f"Probs >> Yes: {prob_yes:.4f}, No: {prob_no:.4f}")
                log("*******************************")
                
                # Update meta-statistics counts and probability sums
                total_fp += fp
                total_fn += fn
                total_tp += tp
                total_tn += tn

                sum_prob_fp += prob_fp
                sum_prob_fn += prob_fn
                sum_prob_tp += prob_tp
                sum_prob_tn += prob_tn

                # Append to true and predicted labels lists for metrics calculation
                true_labels.append(assistant_content)
                predicted_labels.append(generated_response)
                probabilities.append(prob_yes)

                result = {
                    "user_content": user_content,
                    "assistant_content": assistant_content,
                    "generated_response": generated_response,
                    "logprob_yes": logprob_yes,
                    "logprob_no": logprob_no,
                    "prob_yes": prob_yes,
                    "prob_no": prob_no,
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

# Convert labels to binary format for metrics calculation
true_labels_binary = [1 if label == 'Yes' else 0 for label in true_labels]

# Optimal threshold determination
thresholds = np.arange(0.1, 1.0, 0.01)
best_threshold = 0.5
best_f1 = 0

for t in thresholds:
    predicted_labels = (np.array(probabilities) >= t).astype(int)
    f1 = f1_score(true_labels_binary, predicted_labels)
    accuracy = accuracy_score(true_labels_binary, predicted_labels)
    precision = precision_score(true_labels_binary, predicted_labels)
    recall = recall_score(true_labels_binary, predicted_labels)
    log(f"Threshold: {t:.2f}, Accuracy: {accuracy * 100:.2f}%, Precision: {precision * 100:.2f}%, Recall: {recall * 100:.2f}%, F1 Score: {f1 * 100:.2f}%")
    if f1 > best_f1:
        best_f1 = f1
        best_threshold = t

best_threshold = 0.4112 # based on parameters survey
log(f"Optimal Threshold: {best_threshold}")

# Apply the best threshold to make final predictions
final_predictions = (np.array(probabilities) >= best_threshold).astype(int)

# Update the predicted_labels with final predictions for metrics calculation
predicted_labels = ['Yes' if pred == 1 else 'No' for pred in final_predictions]

# Calculate meta-statistics
total_samples = total_fp + total_fn + total_tp + total_tn
meta_statistics = {
    "total_samples": total_samples,
    "total_fp": total_fp,
    "total_fn": total_fn,
    "total_tp": total_tp,
    "total_tn": total_tn,
    "avg_prob_fp": sum_prob_fp / total_fp if total_fp > 0 else 0,
    "avg_prob_fn": sum_prob_fn / total_fn if total_fn > 0 else 0,
    "avg_prob_tp": sum_prob_tp / total_tp if total_tp > 0 else 0,
    "avg_prob_tn": sum_prob_tn / total_tn if total_tn > 0 else 0
}

log("\nMeta-Statistics:")
log(json.dumps(meta_statistics, indent=4))

# Calculate accuracy, precision, recall, specificity, and F1 score
accuracy = accuracy_score(true_labels_binary, final_predictions)
precision = precision_score(true_labels_binary, final_predictions)
recall = recall_score(true_labels_binary, final_predictions)
f1 = f1_score(true_labels_binary, final_predictions)

# Specificity calculation
tn, fp, fn, tp = confusion_matrix(true_labels_binary, final_predictions).ravel()
specificity = tn / (tn + fp)

metrics = {
    "Accuracy": accuracy * 100,
    "Precision": precision * 100,
    "Recall": recall * 100,
    "Specificity": specificity * 100,
    "F1 Score": f1 * 100
}

log("\nMetrics:")
log(json.dumps(metrics, indent=4))

# Save meta-statistics and metrics in a separate file
meta_stats_file = os.path.splitext(output_file_path)[0] + "_meta_stats.json"
with open(meta_stats_file, 'w') as metafile:
    json.dump({"meta_statistics": meta_statistics, "metrics": metrics}, metafile)

log(f"Meta-statistics and metrics saved in {meta_stats_file}")

print(f"Accuracy: {metrics['Accuracy']:.2f}%")
print(f"Precision: {metrics['Precision']:.2f}%")
print(f"Recall: {metrics['Recall']:.2f}%")
print(f"Specificity: {metrics['Specificity']:.2f}%")
print(f"F1 Score: {metrics['F1 Score']:.2f}%")

