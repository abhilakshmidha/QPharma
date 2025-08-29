import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import base64
from io import BytesIO
from typing import Dict, List

def create_summary_plots(processing_times: Dict, file_names: List[str]) -> Dict[str, str]:
    """
    Create summary plots for processing times and return as base64 encoded images
    
    Args:
        processing_times: Dictionary containing 'classical' and/or 'quantum' timing data
        file_names: List of file names that were processed
        
    Returns:
        Dictionary with base64 encoded plot images
    """
    plots = {}
    
    try:
        # Extract timing data
        classical_times = processing_times.get('classical', {})
        quantum_times = processing_times.get('quantum', {})
        
        if not classical_times and not quantum_times:
            return {"error": "No timing data available"}
        
        # Create figure with subplots
        fig_width = 18
        fig_height = 6
        
        # Determine number of plots based on available data
        num_plots = 0
        if classical_times and quantum_times:
            num_plots = 3  # pie chart, bar chart, line graph
        elif classical_times or quantum_times:
            num_plots = 1  # only bar chart for single method
        
        if num_plots == 0:
            return {"error": "No valid timing data"}
        
        fig, axes = plt.subplots(1, num_plots, figsize=(fig_width, fig_height))
        if num_plots == 1:
            axes = [axes]  # Make it iterable
        
        plot_index = 0
        
        # If we have both classical and quantum data
        if classical_times and quantum_times:
            # ==== Pie Chart (Total Execution Time Distribution) ====
            total_classical = sum(classical_times.values())
            total_quantum = sum(quantum_times.values())
            
            methods = ["Classical", "Quantum"]
            times = [total_classical, total_quantum]
            
            axes[plot_index].pie(times, labels=methods, autopct='%1.1f%%', 
                               colors=['blue', 'green'], startangle=140, 
                               explode=(0.1, 0))
            axes[plot_index].set_title("Total Execution Time Distribution")
            plot_index += 1
            
            # ==== Bar Chart (Execution Time Per File) ====
            x = np.arange(len(file_names))
            width = 0.35
            
            classical_values = [classical_times.get(f, 0) for f in file_names]
            quantum_values = [quantum_times.get(f, 0) for f in file_names]
            
            axes[plot_index].bar(x - width/2, classical_values, width, 
                               label="Classical", color='blue', alpha=0.8)
            axes[plot_index].bar(x + width/2, quantum_values, width, 
                               label="Quantum", color='green', alpha=0.8)
            
            axes[plot_index].set_xlabel('Files')
            axes[plot_index].set_ylabel('Time (seconds)')
            axes[plot_index].set_title('Execution Time Comparison per File')
            axes[plot_index].set_xticks(x)
            axes[plot_index].set_xticklabels([f[:10] + '...' if len(f) > 10 else f 
                                           for f in file_names], rotation=45, ha='right')
            axes[plot_index].legend()
            axes[plot_index].grid(True, alpha=0.3)
            plot_index += 1
            
            # ==== Line Graph (Execution Time Trend Across Files) ====
            axes[plot_index].plot(range(len(file_names)), classical_values, 
                                marker='o', linestyle='-', color='blue', 
                                label="Classical", linewidth=2)
            axes[plot_index].plot(range(len(file_names)), quantum_values, 
                                marker='s', linestyle='--', color='green', 
                                label="Quantum", linewidth=2)
            
            axes[plot_index].set_xlabel('File Index')
            axes[plot_index].set_ylabel('Time (seconds)')
            axes[plot_index].set_title('Execution Time Trend Across Files')
            axes[plot_index].legend()
            axes[plot_index].grid(True, alpha=0.3)
            
        else:
            # Single method - just show bar chart
            times_data = classical_times if classical_times else quantum_times
            method_name = "Classical" if classical_times else "Quantum"
            color = 'blue' if classical_times else 'green'
            
            values = [times_data.get(f, 0) for f in file_names]
            
            axes[0].bar(range(len(file_names)), values, color=color, alpha=0.8)
            axes[0].set_xlabel('Files')
            axes[0].set_ylabel('Time (seconds)')
            axes[0].set_title(f'{method_name} Execution Time per File')
            axes[0].set_xticks(range(len(file_names)))
            axes[0].set_xticklabels([f[:10] + '...' if len(f) > 10 else f 
                                   for f in file_names], rotation=45, ha='right')
            axes[0].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Convert plot to base64
        buffer = BytesIO()
        plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
        buffer.seek(0)
        plot_b64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
        buffer.close()
        plt.close()
        
        plots['summary_plot'] = plot_b64
        
        # Generate statistics
        stats = {}
        if classical_times:
            stats['classical'] = {
                'total_time': sum(classical_times.values()),
                'average_time': np.mean(list(classical_times.values())),
                'min_time': min(classical_times.values()),
                'max_time': max(classical_times.values())
            }
        
        if quantum_times:
            stats['quantum'] = {
                'total_time': sum(quantum_times.values()),
                'average_time': np.mean(list(quantum_times.values())),
                'min_time': min(quantum_times.values()),
                'max_time': max(quantum_times.values())
            }
        
        plots['statistics'] = stats
        
    except Exception as e:
        plots['error'] = f"Error creating plots: {str(e)}"
    
    return plots
