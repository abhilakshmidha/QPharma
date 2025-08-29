from fastapi import FastAPI, File, UploadFile, HTTPException, Form
from fastapi.responses import JSONResponse, FileResponse
from fastapi.middleware.cors import CORSMiddleware
import tempfile
import os
import time
import shutil
from typing import List, Optional
import json
import base64
from io import BytesIO
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# Import our processing functions
from classical_processor_demo import process_classical
from quantum_processor_demo import process_quantum
from visualization_utils_demo import create_summary_plots

app = FastAPI(title="Molecular Simulation API", version="1.0.0")

# Add CORS middleware to allow frontend access
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure this properly for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Storage for results (in production, use a database)
results_storage = {}

@app.get("/")
async def root():
    """Root endpoint with API information"""
    return {
        "message": "Molecular Simulation API",
        "endpoints": {
            "upload": "/upload-and-process",
            "classical": "/process-classical",
            "quantum": "/process-quantum",
            "results": "/results/{session_id}",
            "summary": "/summary/{session_id}"
        }
    }

@app.post("/upload-and-process")
async def upload_and_process(
    files: List[UploadFile] = File(...),
    approach: str = Form(..., description="Choose 'classical', 'quantum', or 'both'")
):
    """
    Upload MOL2 files and process them with the specified approach
    """
    if approach not in ["classical", "quantum", "both"]:
        raise HTTPException(
            status_code=400, 
            detail="Approach must be 'classical', 'quantum', or 'both'"
        )
    
    # Validate file types
    for file in files:
        if not file.filename.endswith('.mol2'):
            raise HTTPException(
                status_code=400,
                detail=f"File {file.filename} is not a MOL2 file"
            )
    
    # Create session ID
    session_id = f"session_{int(time.time())}"
    
    # Create temporary directory for this session
    temp_dir = tempfile.mkdtemp(prefix=f"{session_id}_")
    
    try:
        # Save uploaded files
        file_paths = []
        for file in files:
            file_path = os.path.join(temp_dir, file.filename)
            with open(file_path, "wb") as buffer:
                shutil.copyfileobj(file.file, buffer)
            file_paths.append(file_path)
        
        results = {
            "session_id": session_id,
            "files_processed": [f.filename for f in files],
            "approach": approach,
            "classical_results": {},
            "quantum_results": {},
            "processing_time": {},
            "temp_dir": temp_dir
        }
        
        # Process files based on approach
        if approach in ["classical", "both"]:
            classical_results, classical_times = process_classical(file_paths)
            results["classical_results"] = classical_results
            results["processing_time"]["classical"] = classical_times
        
        if approach in ["quantum", "both"]:
            quantum_results, quantum_times = process_quantum(file_paths)
            results["quantum_results"] = quantum_results
            results["processing_time"]["quantum"] = quantum_times
        
        # Store results
        results_storage[session_id] = results
        
        return JSONResponse(content={
            "session_id": session_id,
            "status": "success",
            "message": f"Successfully processed {len(files)} files using {approach} approach(s)",
            "results": {
                "classical": results["classical_results"] if "classical_results" in results else {},
                "quantum": results["quantum_results"] if "quantum_results" in results else {},
                "processing_times": results["processing_time"]
            }
        })
    
    except Exception as e:
        # Clean up temp directory on error
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/process-classical")
async def process_classical_endpoint(files: List[UploadFile] = File(...)):
    """Process files using classical approach only"""
    return await upload_and_process(files, "classical")

@app.post("/process-quantum")
async def process_quantum_endpoint(files: List[UploadFile] = File(...)):
    """Process files using quantum approach only"""
    return await upload_and_process(files, "quantum")

@app.get("/results/{session_id}")
async def get_results(session_id: str):
    """Get detailed results for a session"""
    if session_id not in results_storage:
        raise HTTPException(status_code=404, detail="Session not found")
    
    results = results_storage[session_id]
    return JSONResponse(content={
        "session_id": session_id,
        "files_processed": results["files_processed"],
        "approach": results["approach"],
        "results": {
            "classical": results.get("classical_results", {}),
            "quantum": results.get("quantum_results", {}),
            "processing_times": results.get("processing_time", {})
        }
    })

@app.get("/summary/{session_id}")
async def get_summary(session_id: str):
    """Generate and return summary plots for a session"""
    if session_id not in results_storage:
        raise HTTPException(status_code=404, detail="Session not found")
    
    results = results_storage[session_id]
    
    try:
        # Create summary plots
        plot_data = create_summary_plots(
            results.get("processing_time", {}),
            results["files_processed"]
        )
        
        return JSONResponse(content={
            "session_id": session_id,
            "summary_plots": plot_data,
            "total_files": len(results["files_processed"]),
            "approach": results["approach"]
        })
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error generating summary: {str(e)}")

@app.delete("/cleanup/{session_id}")
async def cleanup_session(session_id: str):
    """Clean up session data and temporary files"""
    if session_id not in results_storage:
        raise HTTPException(status_code=404, detail="Session not found")
    
    results = results_storage[session_id]
    temp_dir = results.get("temp_dir")
    
    # Clean up temporary directory
    if temp_dir and os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    
    # Remove from storage
    del results_storage[session_id]
    
    return JSONResponse(content={
        "message": f"Session {session_id} cleaned up successfully"
    })

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return {"status": "healthy", "timestamp": time.time()}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
