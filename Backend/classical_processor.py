import time
import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
import base64
from typing import List, Dict, Tuple

def visualize_molecule_html(mol, title: str) -> str:
    """Generate HTML visualization of a molecule and return as base64 encoded string"""
    try:
        mol_block = Chem.MolToMolBlock(mol)
        viewer = py3Dmol.view(width=400, height=400)
        viewer.addModel(mol_block, "mol")
        viewer.setStyle({"stick": {}})
        viewer.setBackgroundColor("white")
        viewer.zoomTo()
        
        # Generate HTML
        html_content = viewer._make_html()
        
        # Convert to base64 for easy transmission
        html_b64 = base64.b64encode(html_content.encode('utf-8')).decode('utf-8')
        
        return html_b64
    except Exception as e:
        print(f"Error generating visualization for {title}: {str(e)}")
        return None

def process_classical(file_paths: List[str]) -> Tuple[Dict, Dict]:
    """
    Process MOL2 files using classical computing approach
    
    Args:
        file_paths: List of file paths to process
        
    Returns:
        Tuple of (results_dict, timing_dict)
    """
    classical_results = {}
    classical_times = {}
    
    for file_path in file_paths:
        filename = file_path.split('\\')[-1]  # Get filename from path
        print(f"Processing classical: {filename}")
        
        try:
            # Load molecule
            mol = Chem.MolFromMol2File(file_path, removeHs=False)
            if mol is None:
                classical_results[filename] = {
                    "error": f"Error loading {filename}",
                    "success": False
                }
                classical_times[filename] = 0
                continue
            
            # Classical Simulation: Start Timing
            start_time = time.time()
            
            # Embed molecule with ETKDG method
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            
            # Generate visualization
            visualization_html = visualize_molecule_html(mol, f"Classical: {filename}")
            
            classical_time = time.time() - start_time
            
            # Store results
            classical_results[filename] = {
                "success": True,
                "processing_time": classical_time,
                "visualization": visualization_html,
                "molecule_info": {
                    "num_atoms": mol.GetNumAtoms(),
                    "num_bonds": mol.GetNumBonds(),
                    "molecular_weight": Chem.rdMolDescriptors.CalcExactMolWt(mol)
                }
            }
            
            classical_times[filename] = classical_time
            
            print(f"✅ {filename} - Classical Time: {classical_time:.4f}s")
            
        except Exception as e:
            classical_results[filename] = {
                "error": str(e),
                "success": False
            }
            classical_times[filename] = 0
            print(f"❌ Error processing {filename}: {str(e)}")
    
    return classical_results, classical_times
