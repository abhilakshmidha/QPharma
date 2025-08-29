import time
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
from upload_mol2_files import upload_mol2_files
import webbrowser
import os
import tempfile

def visualize_molecule(mol, title):
    """Render 3D visualization of a molecule"""
    mol_block = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(mol_block, "mol")
    viewer.setStyle({"stick": {}})
    viewer.setBackgroundColor("white")
    viewer.zoomTo()
    print(f"3D Visualization: {title}")
    
    # Save visualization to a temporary HTML file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as temp_file:
        temp_file.write(viewer._make_html().encode("utf-8"))
        temp_file_path = temp_file.name
    
    # Open the HTML file in the default web browser
    webbrowser.open(f"file://{temp_file_path}")


mol2_filenames = upload_mol2_files()

# Process Each File for Classical Computing
classical_times = {}
for mol2_filename in mol2_filenames:
    print(f"\nProcessing: {mol2_filename}")

    # Load molecule
    mol = Chem.MolFromMol2File(mol2_filename, removeHs=False)
    if mol is None:
        print(f"❌ Error loading {mol2_filename}")
        continue

    # Classical Simulation: Start Timing
    start_time_classical = time.time()
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    visualize_molecule(mol, f"Classical: {mol2_filename}")
    classical_time = time.time() - start_time_classical

    # Store Results
    classical_times[mol2_filename] = classical_time

    print(f"✅ {mol2_filename} - Classical Time: {classical_time:.4f}s")

# Print Summary
print("\n📌 Classical Method Summary:")
for filename, time_taken in classical_times.items():
    print(f"⏳ {filename}: {time_taken:.4f} seconds")
