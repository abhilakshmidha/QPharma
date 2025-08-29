import time
import py3Dmol
import pennylane as qml
from rdkit import Chem
import webbrowser
import os
import tempfile
from upload_mol2_files import upload_mol2_files

# Quantum Device Setup
dev = qml.device("default.qubit", wires=4)

@qml.qnode(dev)
def quantum_simulation():
    """Simulate quantum operations"""
    qml.Hadamard(wires=0)
    qml.CNOT(wires=[0, 1])
    qml.RX(0.5, wires=2)
    qml.CNOT(wires=[2, 3])
    return qml.probs(wires=[0, 1, 2, 3])

mol2_filenames = upload_mol2_files()

# Visualization Function
def visualize_molecule(mol, title):
    """Render 3D visualization of a molecule"""
    mol_block = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({"stick": {}})
    viewer.setBackgroundColor("white")
    viewer.zoomTo()
    print(f"3D Visualization: {title}")
    # return viewer.show()
    # Save visualization to a temporary HTML file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as temp_file:
        temp_file.write(viewer._make_html().encode("utf-8"))
        temp_file_path = temp_file.name
    
    # Open the HTML file in the default web browser
    webbrowser.open(f"file://{temp_file_path}")

quantum_times = {}
quantum_results = {}

for mol2_filename in mol2_filenames:
    print(f"\nProcessing: {mol2_filename}")

    # Load molecule
    mol = Chem.MolFromMol2File(mol2_filename, removeHs=False)
    if mol is None:
        print(f"❌ Error loading {mol2_filename}")
        continue

    # Quantum Simulation: Start Timing
    start_time_quantum = time.time()
    quantum_result = quantum_simulation()
    visualize_molecule(mol, f"Quantum: {mol2_filename}")
    quantum_time = time.time() - start_time_quantum

    # Store Results
    quantum_times[mol2_filename] = quantum_time
    quantum_results[mol2_filename] = quantum_result

    print(f"✅ {mol2_filename} - Quantum Time: {quantum_time:.4f}s")
    print(f"Quantum Simulation Output: {quantum_result}")

# Print Summary
print("\n📌 Quantum Method Summary:")
for filename, time_taken in quantum_times.items():
    print(f"⏳ {filename}: {time_taken:.4f} seconds")