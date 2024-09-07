import sys
import matplotlib

# Use the 'Agg' backend to avoid tkinter dependency issues
matplotlib.use('Agg')

import pandas as pd
import matplotlib.pyplot as plt
import os

def plot_trends(filepaths):
    # Define colors for different files
    colors = ['blue', 'orange', 'green']
    
    # Initialize figures
    plt.figure(figsize=(12, 6))
    
    # Plot Train and Valid Loss
    plt.subplot(1, 2, 1)
    for i, filepath in enumerate(filepaths):
        data = pd.read_csv(filepath)
        data = data[data['step'] <= 914]  # Filter steps up to 914
        filename = os.path.splitext(os.path.basename(filepath))[0]
        plt.plot(data['step'], data['train_loss'], label=f'Train Loss {filename}', color=colors[i], linestyle='-', linewidth=0.5)
        plt.plot(data['step'], data['valid_loss'], label=f'Valid Loss {filename}', color=colors[i], linestyle='--', linewidth=0.5)
    plt.xlabel('Step')
    plt.ylabel('Loss')
    plt.title('Train and Valid Loss')
    plt.legend()
    
    # Plot Train and Valid Mean Token Accuracy
    plt.subplot(1, 2, 2)
    for i, filepath in enumerate(filepaths):
        data = pd.read_csv(filepath)
        data = data[data['step'] <= 914]  # Filter steps up to 914
        filename = os.path.splitext(os.path.basename(filepath))[0]
        plt.plot(data['step'], data['train_mean_token_accuracy'], label=f'Train Accuracy {filename}', color=colors[i], linestyle='-', linewidth=0.5)
        plt.plot(data['step'], data['valid_mean_token_accuracy'], label=f'Valid Accuracy {filename}', color=colors[i], linestyle='--', linewidth=0.5)
    plt.xlabel('Step')
    plt.ylabel('Mean Token Accuracy')
    plt.title('Train and Valid Mean Token Accuracy')
    plt.legend()

    plt.tight_layout()
    plt.savefig('train_valid_trends.png')  # Save the figure as a file

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python plot_trends.py <file1.csv> <file2.csv> <file3.csv>")
    else:
        filepaths = sys.argv[1:]
        plot_trends(filepaths)

