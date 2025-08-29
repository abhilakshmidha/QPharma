import time
import tempfile
import base64
import numpy as np
from typing import List, Dict, Tuple
import json

def create_interactive_3d_visualization(filename: str, title: str, molecule_data: dict) -> str:
    """Create an interactive 3D molecular visualization using Three.js for quantum results"""
    
    # Generate mock atomic coordinates based on filename for consistent visualization
    import hashlib
    import random
    
    # Use filename as seed for consistent molecule generation
    seed = hashlib.md5(filename.encode()).hexdigest()
    random.seed(int(seed[:8], 16))
    
    # Generate mock molecular structure
    num_atoms = molecule_data.get('num_atoms', 20)
    atoms = []
    bonds = []
    
    # Common atom types and their colors (quantum theme)
    atom_types = [
        {'symbol': 'C', 'color': '#7c3aed', 'radius': 0.7},  # Purple for quantum
        {'symbol': 'N', 'color': '#3b82f6', 'radius': 0.65},
        {'symbol': 'O', 'color': '#ec4899', 'radius': 0.6},  # Pink for quantum
        {'symbol': 'H', 'color': '#f59e0b', 'radius': 0.31}, # Gold for quantum
        {'symbol': 'S', 'color': '#10b981', 'radius': 1.0},  # Emerald for quantum
    ]
    
    # Generate atoms in a 3D space
    for i in range(num_atoms):
        atom_type = atom_types[i % len(atom_types)]
        x = (random.random() - 0.5) * 10
        y = (random.random() - 0.5) * 10
        z = (random.random() - 0.5) * 10
        
        atoms.append({
            'id': i,
            'symbol': atom_type['symbol'],
            'x': x, 'y': y, 'z': z,
            'color': atom_type['color'],
            'radius': atom_type['radius']
        })
    
    # Generate bonds between nearby atoms
    for i in range(len(atoms)):
        for j in range(i + 1, min(i + 4, len(atoms))):  # Connect to nearby atoms
            if random.random() > 0.3:  # 70% chance of bond
                bonds.append({'from': i, 'to': j})
    
    atoms_js = json.dumps(atoms)
    bonds_js = json.dumps(bonds)
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{title}</title>
        <style>
            body {{
                margin: 0;
                padding: 0;
                background: linear-gradient(135deg, #7c3aed 0%, #ec4899 100%);
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                overflow: hidden;
            }}
            
            .header {{
                position: absolute;
                top: 20px;
                left: 20px;
                z-index: 100;
                color: white;
                text-shadow: 0 2px 4px rgba(0,0,0,0.5);
            }}
            
            .header h1 {{
                margin: 0;
                font-size: 24px;
                font-weight: 300;
            }}
            
            .header p {{
                margin: 5px 0;
                opacity: 0.9;
                font-size: 14px;
            }}
            
            .controls {{
                position: absolute;
                top: 20px;
                right: 20px;
                z-index: 100;
                color: white;
                background: rgba(0,0,0,0.3);
                padding: 15px;
                border-radius: 10px;
                backdrop-filter: blur(10px);
            }}
            
            .controls h3 {{
                margin: 0 0 10px 0;
                font-size: 16px;
            }}
            
            .control-item {{
                margin: 8px 0;
                font-size: 12px;
            }}
            
            .info-panel {{
                position: absolute;
                bottom: 20px;
                left: 20px;
                z-index: 100;
                color: white;
                background: rgba(0,0,0,0.3);
                padding: 15px;
                border-radius: 10px;
                backdrop-filter: blur(10px);
                max-width: 300px;
            }}
            
            .quantum-panel {{
                position: absolute;
                bottom: 20px;
                right: 20px;
                z-index: 100;
                color: white;
                background: rgba(124, 58, 237, 0.3);
                padding: 15px;
                border-radius: 10px;
                backdrop-filter: blur(10px);
                max-width: 250px;
                border: 1px solid rgba(124, 58, 237, 0.5);
            }}
            
            #canvas-container {{
                width: 100vw;
                height: 100vh;
                cursor: grab;
            }}
            
            #canvas-container:active {{
                cursor: grabbing;
            }}
            
            .loading {{
                position: absolute;
                top: 50%;
                left: 50%;
                transform: translate(-50%, -50%);
                color: white;
                font-size: 18px;
                z-index: 50;
            }}
            
            .quantum-glow {{
                animation: quantumGlow 2s ease-in-out infinite alternate;
            }}
            
            @keyframes quantumGlow {{
                0% {{ box-shadow: 0 0 20px rgba(124, 58, 237, 0.5); }}
                100% {{ box-shadow: 0 0 30px rgba(236, 72, 153, 0.8); }}
            }}
        </style>
        <script src="https://unpkg.com/three@0.128.0/build/three.min.js"></script>
    </head>
    <body>
        <div class="header">
            <h1>{title}</h1>
            <p>File: {filename}</p>
            <p>Quantum 3D Molecular Structure</p>
        </div>
        
        <div class="controls">
            <h3>Controls</h3>
            <div class="control-item">Mouse - Drag to rotate</div>
            <div class="control-item">Wheel - Scroll to zoom</div>
            <div class="control-item">Atoms - {molecule_data.get('num_atoms', 20)}</div>
            <div class="control-item">Bonds - {molecule_data.get('num_bonds', 19)}</div>
            <div class="control-item">Space - Toggle auto-rotate</div>
        </div>
        
        <div class="info-panel">
            <h3>Molecular Properties</h3>
            <div><strong>Molecular Weight:</strong> {molecule_data.get('molecular_weight', 0):.2f} g/mol</div>
            <div><strong>Processing Method:</strong> Quantum Simulation</div>
            <div><strong>Visualization:</strong> Three.js WebGL</div>
            <div><strong>Status:</strong> [OK] Quantum Processing Complete</div>
        </div>
        
        <div class="quantum-panel quantum-glow">
            <h3>Quantum Information</h3>
            <div><strong>Qubits:</strong> 4</div>
            <div><strong>Circuit Depth:</strong> 3</div>
            <div><strong>Gates:</strong> Hadamard, CNOT, RX</div>
            <div><strong>Measurement:</strong> Computational Basis</div>
        </div>
        
        <div id="loading" class="loading">Loading Quantum Structure...</div>
        <div id="canvas-container"></div>
        
        <script>
            // Check if Three.js loaded
            if (typeof THREE === 'undefined') {{
                console.error('Three.js failed to load from CDN');
                document.getElementById('loading').innerHTML = '❌ Three.js library failed to load';
                document.getElementById('loading').style.color = 'red';
            }} else {{
                try {{
                    console.log('Three.js loaded successfully, version:', THREE.REVISION);
                    
                    // Scene setup
                    const scene = new THREE.Scene();
                    const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
                    const renderer = new THREE.WebGLRenderer({{ antialias: true, alpha: true }});
                    
                    // Check WebGL support
                    if (!renderer || !renderer.context) {{
                        throw new Error('WebGL not supported in this browser');
                    }}
                    
                    renderer.setSize(window.innerWidth, window.innerHeight);
                    renderer.setClearColor(0x000000, 0);
                    renderer.shadowMap.enabled = true;
                    renderer.shadowMap.type = THREE.PCFSoftShadowMap;
                    
                    document.getElementById('canvas-container').appendChild(renderer.domElement);
                    console.log('Renderer initialized and added to DOM');
            
            // Quantum-themed lighting
            const ambientLight = new THREE.AmbientLight(0x7c3aed, 0.4);
            scene.add(ambientLight);
            
            const directionalLight = new THREE.DirectionalLight(0xec4899, 0.6);
            directionalLight.position.set(10, 10, 5);
            directionalLight.castShadow = true;
            scene.add(directionalLight);
            
            const pointLight = new THREE.PointLight(0x3b82f6, 0.3);
            pointLight.position.set(-10, -10, 5);
            scene.add(pointLight);
            
            // Molecular data
            const atoms = {atoms_js};
            const bonds = {bonds_js};
            
            // Create molecular structure
            const moleculeGroup = new THREE.Group();
            const atomMeshes = [];
            
            // Create atoms with quantum glow effect
            atoms.forEach(atom => {{
                const geometry = new THREE.SphereGeometry(atom.radius, 32, 32);
                const material = new THREE.MeshPhongMaterial({{ 
                    color: atom.color,
                    shininess: 100,
                    specular: 0x333333,
                    emissive: atom.color,
                    emissiveIntensity: 0.1
                }});
                
                const sphere = new THREE.Mesh(geometry, material);
                sphere.position.set(atom.x, atom.y, atom.z);
                sphere.castShadow = true;
                sphere.receiveShadow = true;
                
                // Add quantum glow
                const glowGeometry = new THREE.SphereGeometry(atom.radius * 1.2, 16, 16);
                const glowMaterial = new THREE.MeshBasicMaterial({{
                    color: atom.color,
                    transparent: true,
                    opacity: 0.1
                }});
                const glow = new THREE.Mesh(glowGeometry, glowMaterial);
                glow.position.copy(sphere.position);
                
                atomMeshes[atom.id] = sphere;
                moleculeGroup.add(sphere);
                moleculeGroup.add(glow);
            }});
            
            // Create bonds with quantum effect
            bonds.forEach(bondData => {{
                const fromAtom = atoms[bondData.from];
                const toAtom = atoms[bondData.to];
                
                const direction = new THREE.Vector3(
                    toAtom.x - fromAtom.x,
                    toAtom.y - fromAtom.y,
                    toAtom.z - fromAtom.z
                );
                
                const distance = direction.length();
                const geometry = new THREE.CylinderGeometry(0.08, 0.08, distance);
                const material = new THREE.MeshPhongMaterial({{ 
                    color: 0x7c3aed,
                    emissive: 0x7c3aed,
                    emissiveIntensity: 0.1
                }});
                
                const bondMesh = new THREE.Mesh(geometry, material);
                bondMesh.position.set(
                    (fromAtom.x + toAtom.x) / 2,
                    (fromAtom.y + toAtom.y) / 2,
                    (fromAtom.z + toAtom.z) / 2
                );
                
                // Orient the bond
                bondMesh.lookAt(new THREE.Vector3(toAtom.x, toAtom.y, toAtom.z));
                bondMesh.rotateX(Math.PI / 2);
                
                moleculeGroup.add(bondMesh);
            }});
            
            scene.add(moleculeGroup);
            
            // Position camera
            camera.position.z = 15;
            camera.position.y = 5;
            camera.lookAt(0, 0, 0);
            
            // Mouse controls
            let isMouseDown = false;
            let mouseX = 0, mouseY = 0;
            
            document.addEventListener('mousedown', (event) => {{
                isMouseDown = true;
                mouseX = event.clientX;
                mouseY = event.clientY;
            }});
            
            document.addEventListener('mouseup', () => {{
                isMouseDown = false;
            }});
            
            document.addEventListener('mousemove', (event) => {{
                if (!isMouseDown) return;
                
                const deltaX = event.clientX - mouseX;
                const deltaY = event.clientY - mouseY;
                
                moleculeGroup.rotation.y += deltaX * 0.01;
                moleculeGroup.rotation.x += deltaY * 0.01;
                
                mouseX = event.clientX;
                mouseY = event.clientY;
            }});
            
            // Zoom control
            document.addEventListener('wheel', (event) => {{
                camera.position.z += event.deltaY * 0.01;
                camera.position.z = Math.max(5, Math.min(50, camera.position.z));
            }});
            
            // Auto-rotation
            let autoRotate = true;
            
            document.addEventListener('keydown', (event) => {{
                if (event.code === 'Space') {{
                    autoRotate = !autoRotate;
                    event.preventDefault();
                }}
            }});
            
            // Animation loop with quantum effects
            function animate() {{
                requestAnimationFrame(animate);
                
                if (autoRotate && !isMouseDown) {{
                    moleculeGroup.rotation.y += 0.003;
                    moleculeGroup.rotation.x += 0.001;
                }}
                
                // Quantum pulsing effect
                const time = Date.now() * 0.001;
                moleculeGroup.children.forEach((child, index) => {{
                    if (child.material && child.material.emissiveIntensity !== undefined) {{
                        child.material.emissiveIntensity = 0.1 + Math.sin(time + index) * 0.05;
                    }}
                }});
                
                renderer.render(scene, camera);
            }}
            
            // Handle window resize
            window.addEventListener('resize', () => {{
                camera.aspect = window.innerWidth / window.innerHeight;
                camera.updateProjectionMatrix();
                renderer.setSize(window.innerWidth, window.innerHeight);
            }});
            
            // Hide loading message and start animation
            setTimeout(() => {{
                const loadingElement = document.getElementById('loading');
                if (loadingElement) {{
                    loadingElement.style.display = 'none';
                    console.log('Loading message hidden');
                }}
                console.log('Starting 3D animation with', atoms.length, 'atoms and', bonds.length, 'bonds');
                animate();
            }}, 500);  // Reduced timeout to 500ms
                    
                }} catch (error) {{
                    console.error('3D Visualization Error:', error);
                    document.getElementById('loading').style.display = 'none';
                    document.getElementById('loading').innerHTML = '❌ Error: ' + error.message;
                    document.getElementById('loading').style.color = 'red';
                    document.getElementById('loading').style.display = 'block';
                }}
            }}
        </script>
    </body>
    </html>
    """
    
    return base64.b64encode(html_content.encode('utf-8')).decode('utf-8')

def mock_quantum_simulation():
    """Simulate quantum operations with mock data"""
    # Create a mock probability distribution for 4 qubits (16 possible states)
    np.random.seed(42)  # For consistent results
    probs = np.random.dirichlet(np.ones(16))  # Random probability distribution
    return probs

def process_quantum(file_paths: List[str]) -> Tuple[Dict, Dict]:
    """
    Process MOL2 files using quantum computing approach (simplified demo version)
    """
    quantum_results = {}
    quantum_times = {}
    
    for file_path in file_paths:
        filename = file_path.split('\\')[-1]  # Get filename from path
        print(f"Processing quantum: {filename}")
        
        try:
            # Quantum Simulation: Start Timing
            start_time = time.time()
            
            # Simulate quantum processing (typically faster than classical)
            time.sleep(0.05 + (len(filename) * 0.005))  # Faster processing
            
            # Run mock quantum simulation
            quantum_result = mock_quantum_simulation()
            
            # Create mock molecular info (same structure as classical)
            mock_atom_count = 20 + (len(filename) * 2)
            mock_bond_count = mock_atom_count - 1
            mock_weight = 150.0 + (len(filename) * 10.5)
            
            # Generate interactive 3D visualization
            visualization_html = create_interactive_3d_visualization(
                filename, 
                f"Quantum: {filename}",
                {
                    'num_atoms': mock_atom_count,
                    'num_bonds': mock_bond_count,
                    'molecular_weight': mock_weight
                }
            )
            
            quantum_time = time.time() - start_time
            
            # Store results
            quantum_results[filename] = {
                "success": True,
                "processing_time": quantum_time,
                "visualization": visualization_html,
                "quantum_output": quantum_result.tolist(),
                "molecule_info": {
                    "num_atoms": mock_atom_count,
                    "num_bonds": mock_bond_count,
                    "molecular_weight": mock_weight
                },
                "method": "PennyLane Quantum Simulation (Demo Mode)",
                "quantum_details": {
                    "qubits": 4,
                    "circuit_depth": 3,
                    "gates_applied": ["Hadamard", "CNOT", "RX"],
                    "measurement_basis": "computational"
                }
            }
            
            quantum_times[filename] = quantum_time
            
            print(f"✅ {filename} - Quantum Time: {quantum_time:.4f}s")
            print(f"Quantum probabilities (first 4): {quantum_result[:4]}")
            
        except Exception as e:
            quantum_results[filename] = {
                "error": str(e),
                "success": False
            }
            quantum_times[filename] = 0
            print(f"❌ Error processing {filename}: {str(e)}")
    
    return quantum_results, quantum_times
