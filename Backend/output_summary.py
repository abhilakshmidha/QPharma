import matplotlib.pyplot as plt
import numpy as np

# Assuming ⁠ classical_times ⁠ and ⁠ quantum_times ⁠ are already defined as lists
# Example values (Replace with actual execution times)
# classical_times = [12, 15, 13, 16, 14]
# quantum_times = [3, 4, 3.5, 4.2, 3.8]

# Number of files (assuming same length for both lists)
num_files = len(classical_times)
files = [f"File {i+1}" for i in range(num_files)]  # Generate file labels dynamically

# Summing total execution time
total_classical_time = sum(classical_times.values())
total_quantum_time = sum(quantum_times.values())

# Labels and values for pie chart
methods = ["Classical", "Quantum"]
times = [total_classical_time, total_quantum_time]

# Create a figure with subplots
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# ==== Pie Chart (Total Execution Time Distribution) ====
axes[0].pie(times, labels=methods, autopct='%1.1f%%', colors=['blue', 'green'], startangle=140, explode=(0.1, 0))
axes[0].set_title("Total Execution Time Distribution")

# ==== Bar Chart (Execution Time Per File) ====
x = np.arange(num_files)  # X-axis positions
width = 0.4  # Bar width

axes[1].bar(x - width/2, classical_times.values(), width, label="Classical", color='blue')
axes[1].bar(x + width/2, quantum_times.values(), width, label="Quantum", color='green')

axes[1].set_xticks(x)
axes[1].set_xticklabels(files, rotation=45, ha="right")
axes[1].set_ylabel("Time (seconds)")
axes[1].set_title("Execution Time Comparison per File")
axes[1].legend()

# ==== Line Graph (Execution Time Trend Across Files) ====
axes[2].plot(files, classical_times.values(), marker='o', linestyle='-', color='blue', label="Classical")
axes[2].plot(files, quantum_times.values(), marker='s', linestyle='--', color='green', label="Quantum")

axes[2].set_ylabel("Time (seconds)")
axes[2].set_title("Execution Time Trend Across Files")
axes[2].legend()

# Show all plots
plt.tight_layout()
plt.show()