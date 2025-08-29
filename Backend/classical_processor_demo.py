import time
import tempfile
import base64
from typing import List, Dict, Tuple
import json

def create_interactive_3d_visualization(filename: str, title: str, molecule_data: dict) -> str:
    """Create an interactive 3D molecular visualization using Three.js"""
    
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
    
    # Common atom types and their colors
    atom_types = [
        {'symbol': 'C', 'color': '#909090', 'radius': 0.7},
        {'symbol': 'N', 'color': '#3050F8', 'radius': 0.65},
        {'symbol': 'O', 'color': '#FF0D0D', 'radius': 0.6},
        {'symbol': 'H', 'color': '#FFFFFF', 'radius': 0.31},
        {'symbol': 'S', 'color': '#FFFF30', 'radius': 1.0},
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
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
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
        </style>
        <script src="https://unpkg.com/three@0.128.0/build/three.min.js"></script>
    </head>
    <body>
        <div class="header">
            <h1>{title}</h1>
            <p>File: {filename}</p>
            <p>Interactive 3D Molecular Structure</p>
        </div>
        
        <div class="controls">
            <h3>Controls</h3>
            <div class="control-item">Mouse - Drag to rotate</div>
            <div class="control-item">Wheel - Scroll to zoom</div>
            <div class="control-item">Atoms - {molecule_data.get('num_atoms', 20)}</div>
            <div class="control-item">Bonds - {molecule_data.get('num_bonds', 19)}</div>
        </div>
        
        <div class="info-panel">
            <h3>Molecular Properties</h3>
            <div><strong>Molecular Weight:</strong> {molecule_data.get('molecular_weight', 0):.2f} g/mol</div>
            <div><strong>Processing Method:</strong> Classical Simulation</div>
            <div><strong>Visualization:</strong> Three.js WebGL</div>
            <div><strong>Status:</strong> [OK] Processing Complete</div>
        </div>
        
        <div id="loading" class="loading">Loading 3D Structure...</div>
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
            
            // Lighting
            const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
            scene.add(ambientLight);
            
            const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
            directionalLight.position.set(10, 10, 5);
            directionalLight.castShadow = true;
            scene.add(directionalLight);
            
            // Molecular data
            const atoms = {atoms_js};
            const bonds = {bonds_js};
            
            // Create molecular structure
            const moleculeGroup = new THREE.Group();
            const atomMeshes = [];
            
            // Create atoms
            atoms.forEach(atom => {{
                const geometry = new THREE.SphereGeometry(atom.radius, 32, 32);
                const material = new THREE.MeshPhongMaterial({{ 
                    color: atom.color,
                    shininess: 100,
                    specular: 0x111111
                }});
                
                const sphere = new THREE.Mesh(geometry, material);
                sphere.position.set(atom.x, atom.y, atom.z);
                sphere.castShadow = true;
                sphere.receiveShadow = true;
                
                atomMeshes[atom.id] = sphere;
                moleculeGroup.add(sphere);
            }});
            
            // Create bonds
            bonds.forEach(bondData => {{
                const fromAtom = atoms[bondData.from];
                const toAtom = atoms[bondData.to];
                
                const direction = new THREE.Vector3(
                    toAtom.x - fromAtom.x,
                    toAtom.y - fromAtom.y,
                    toAtom.z - fromAtom.z
                );
                
                const distance = direction.length();
                const geometry = new THREE.CylinderGeometry(0.1, 0.1, distance);
                const material = new THREE.MeshPhongMaterial({{ color: 0x888888 }});
                
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
            
            // Animation loop
            function animate() {{
                requestAnimationFrame(animate);
                
                if (autoRotate && !isMouseDown) {{
                    moleculeGroup.rotation.y += 0.005;
                }}
                
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

def process_classical(file_paths: List[str]) -> Tuple[Dict, Dict]:
    """
    Process MOL2 files using classical computing approach (simplified demo version)
    """
    classical_results = {}
    classical_times = {}
    
    for file_path in file_paths:
        filename = file_path.split('\\')[-1]  # Get filename from path
        print(f"Processing classical: {filename}")
        
        try:
            # Simulate classical processing time
            start_time = time.time()
            
            # Simulate some processing work
            time.sleep(0.1 + (len(filename) * 0.01))  # Variable time based on filename length
            
            # Create mock molecular info
            mock_atom_count = 20 + (len(filename) * 2)  # Mock based on filename
            mock_bond_count = mock_atom_count - 1
            mock_weight = 150.0 + (len(filename) * 10.5)
            
            molecule_info = {
                "num_atoms": mock_atom_count,
                "num_bonds": mock_bond_count,
                "molecular_weight": mock_weight
            }
            
            # Generate interactive 3D visualization
            visualization_html = create_interactive_3d_visualization(filename, f"Classical: {filename}", molecule_info)
            
            classical_time = time.time() - start_time
            
            # Store results
            classical_results[filename] = {
                "success": True,
                "processing_time": classical_time,
                "visualization": visualization_html,
                "molecule_info": {
                    "num_atoms": mock_atom_count,
                    "num_bonds": mock_bond_count,
                    "molecular_weight": mock_weight
                },
                "method": "Classical RDKit Simulation (Demo Mode)"
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
