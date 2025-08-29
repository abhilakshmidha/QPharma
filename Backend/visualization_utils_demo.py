import base64
from io import BytesIO
from typing import Dict, List
import json
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np

def create_summary_plots(processing_times: Dict, file_names: List[str]) -> Dict[str, str]:
    """
    Create summary plots for processing times using matplotlib
    """
    plots = {}
    
    try:
        # Extract timing data
        classical_times = processing_times.get('classical', {})
        quantum_times = processing_times.get('quantum', {})
        
        if not classical_times and not quantum_times:
            return {"error": "No timing data available"}
        
        # Create matplotlib figure
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        fig.patch.set_facecolor('white')
        
        # Get all files that were processed
        all_files = list(set(list(classical_times.keys()) + list(quantum_times.keys())))
        
        # Prepare data for charts
        classical_values = [classical_times.get(f, 0) for f in all_files]
        quantum_values = [quantum_times.get(f, 0) for f in all_files]
        
        # 1. Pie Chart - Total execution time distribution
        if classical_times or quantum_times:
            total_classical = sum(classical_times.values()) if classical_times else 0
            total_quantum = sum(quantum_times.values()) if quantum_times else 0
            
            if total_classical > 0 and total_quantum > 0:
                methods = ['Classical', 'Quantum']
                totals = [total_classical, total_quantum]
                colors = ['#3b82f6', '#7c3aed']
                
                wedges, texts, autotexts = axes[0].pie(totals, labels=methods, autopct='%1.1f%%', 
                                                      colors=colors, startangle=140, explode=(0.1, 0))
                axes[0].set_title('Total Execution Time Distribution', fontsize=14, fontweight='bold')
                
                # Style the text
                for autotext in autotexts:
                    autotext.set_color('white')
                    autotext.set_fontweight('bold')
            else:
                axes[0].text(0.5, 0.5, 'Insufficient data\\nfor pie chart', 
                           ha='center', va='center', transform=axes[0].transAxes, fontsize=12)
                axes[0].set_title('Total Execution Time Distribution', fontsize=14, fontweight='bold')
        
        # 2. Bar Chart - Execution time per file
        if all_files:
            x = np.arange(len(all_files))
            width = 0.35
            
            bars1 = bars2 = None
            if classical_times:
                bars1 = axes[1].bar(x - width/2, classical_values, width, 
                                   label='Classical', color='#3b82f6', alpha=0.8)
            if quantum_times:
                bars2 = axes[1].bar(x + width/2, quantum_values, width, 
                                   label='Quantum', color='#7c3aed', alpha=0.8)
            
            axes[1].set_xlabel('Files')
            axes[1].set_ylabel('Time (seconds)')
            axes[1].set_title('Execution Time Comparison per File', fontsize=14, fontweight='bold')
            axes[1].set_xticks(x)
            axes[1].set_xticklabels([f[:8] + '...' if len(f) > 8 else f for f in all_files], 
                                   rotation=45, ha='right')
            axes[1].legend()
            axes[1].grid(axis='y', alpha=0.3)
            
            # Add value labels on bars
            if bars1:
                for bar, val in zip(bars1, classical_values):
                    if val > 0:
                        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001, 
                                   f'{val:.3f}s', ha='center', va='bottom', fontsize=8)
            if bars2:
                for bar, val in zip(bars2, quantum_values):
                    if val > 0:
                        axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001, 
                                   f'{val:.3f}s', ha='center', va='bottom', fontsize=8)
        
        # 3. Line Graph - Execution time trend
        if all_files:
            x_pos = range(len(all_files))
            
            if classical_times:
                axes[2].plot(x_pos, classical_values, marker='o', linestyle='-', 
                           color='#3b82f6', label='Classical', linewidth=2, markersize=6)
            if quantum_times:
                axes[2].plot(x_pos, quantum_values, marker='s', linestyle='--', 
                           color='#7c3aed', label='Quantum', linewidth=2, markersize=6)
            
            axes[2].set_xlabel('Files')
            axes[2].set_ylabel('Time (seconds)')
            axes[2].set_title('Execution Time Trend Across Files', fontsize=14, fontweight='bold')
            axes[2].set_xticks(x_pos)
            axes[2].set_xticklabels([f[:8] + '...' if len(f) > 8 else f for f in all_files], 
                                   rotation=45, ha='right')
            axes[2].legend()
            axes[2].grid(alpha=0.3)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save to base64
        buffer = BytesIO()
        plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight', facecolor='white')
        buffer.seek(0)
        plot_b64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
        plt.close(fig)
        
        plots['summary_plot'] = plot_b64
        
        # Generate statistics
        stats = {}
        if classical_times:
            times_list = list(classical_times.values())
            stats['classical'] = {
                'total_time': sum(times_list),
                'average_time': sum(times_list) / len(times_list),
                'min_time': min(times_list),
                'max_time': max(times_list)
            }
        
        if quantum_times:
            times_list = list(quantum_times.values())
            stats['quantum'] = {
                'total_time': sum(times_list),
                'average_time': sum(times_list) / len(times_list),
                'min_time': min(times_list),
                'max_time': max(times_list)
            }
        
        plots['statistics'] = stats
        
    except Exception as e:
        plots['error'] = f"Error creating plots: {str(e)}"
        print(f"Error in create_summary_plots: {e}")
    
    return plots

