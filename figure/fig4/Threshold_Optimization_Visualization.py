import matplotlib.pyplot as plt
import sys

# Function to parse the data from the text file
def parse_data(file_path):
    thresholds = []
    accuracies = []
    precisions = []
    recalls = []
    f1_scores = []

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            parts = line.strip().split(',')
            threshold = float(parts[0].split(': ')[1])
            accuracy = float(parts[1].split(': ')[1].strip('%')) / 100.0
            precision = float(parts[2].split(': ')[1].strip('%')) / 100.0
            recall = float(parts[3].split(': ')[1].strip('%')) / 100.0
            f1_score = float(parts[4].split(': ')[1].strip('%')) / 100.0

            thresholds.append(threshold)
            accuracies.append(accuracy)
            precisions.append(precision)
            recalls.append(recall)
            f1_scores.append(f1_score)

    return thresholds, accuracies, precisions, recalls, f1_scores

# File path to your data
input_file_path = sys.argv[1]

# Parse data from the file
thresholds, accuracies, precisions, recalls, f1_scores = parse_data(file_path)

# Plotting
plt.figure(figsize=(10, 6))

plt.plot(thresholds, accuracies, label='Accuracy', marker='o')
plt.plot(thresholds, precisions, label='Precision', marker='o')
plt.plot(thresholds, recalls, label='Recall', marker='o')
plt.plot(thresholds, f1_scores, label='F1 Score', marker='o')

plt.xlabel('Threshold')
plt.ylabel('Metrics')
plt.title('Metrics vs. Threshold')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Show plot
plt.show()

