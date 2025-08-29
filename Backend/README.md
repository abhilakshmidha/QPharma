# Molecular Simulation FastAPI Backend

A FastAPI-based backend for processing MOL2 molecular files using both classical and quantum computing approaches.

## Features

- **File Upload**: Upload multiple MOL2 files via REST API
- **Classical Processing**: Process molecules using RDKit and classical computing methods
- **Quantum Processing**: Process molecules using PennyLane quantum computing framework
- **3D Visualization**: Generate interactive 3D molecular visualizations using py3Dmol
- **Performance Analysis**: Compare processing times between classical and quantum approaches
- **Summary Reports**: Generate statistical summaries and visualizations

## API Endpoints

### Core Endpoints

- `POST /upload-and-process` - Upload files and process with specified approach
- `POST /process-classical` - Process files using classical approach only
- `POST /process-quantum` - Process files using quantum approach only
- `GET /results/{session_id}` - Get detailed results for a processing session
- `GET /summary/{session_id}` - Get summary plots and statistics
- `DELETE /cleanup/{session_id}` - Clean up session data and temporary files
- `GET /health` - Health check endpoint

## Installation

1. **Install Python dependencies:**

   ```cmd
   pip install -r requirements.txt
   ```

2. **Start the server:**

   ```cmd
   # Option 1: Use the provided batch script
   start_server.bat

   # Option 2: Start manually
   uvicorn main:app --host 0.0.0.0 --port 8000 --reload
   ```

## Usage

### 1. Start the Server

Run `start_server.bat` or use uvicorn directly. The server will be available at `http://localhost:8000`.

### 2. Upload and Process Files

**Using curl:**

```bash
curl -X POST "http://localhost:8000/upload-and-process" \
  -H "Content-Type: multipart/form-data" \
  -F "files=@117_ideal.mol2" \
  -F "approach=both"
```

**Using Python requests:**

```python
import requests

files = [('files', ('117_ideal.mol2', open('117_ideal.mol2', 'rb')))]
data = {'approach': 'both'}  # 'classical', 'quantum', or 'both'

response = requests.post('http://localhost:8000/upload-and-process',
                        files=files, data=data)
print(response.json())
```

### 3. Get Results

```python
session_id = "your_session_id_here"
response = requests.get(f'http://localhost:8000/results/{session_id}')
print(response.json())
```

### 4. Get Summary Plots

```python
response = requests.get(f'http://localhost:8000/summary/{session_id}')
result = response.json()

# Save the plot
import base64
with open('summary_plot.png', 'wb') as f:
    f.write(base64.b64decode(result['summary_plots']['summary_plot']))
```

## API Response Format

### Upload Response

```json
{
  "session_id": "session_1693123456",
  "status": "success",
  "message": "Successfully processed 1 files using both approach(s)",
  "results": {
    "classical": {
      "117_ideal.mol2": {
        "success": true,
        "processing_time": 0.1234,
        "visualization": "base64_encoded_html",
        "molecule_info": {
          "num_atoms": 25,
          "num_bonds": 24,
          "molecular_weight": 180.157
        }
      }
    },
    "quantum": {
      "117_ideal.mol2": {
        "success": true,
        "processing_time": 0.0567,
        "visualization": "base64_encoded_html",
        "quantum_output": [0.25, 0.25, 0.25, 0.25],
        "molecule_info": {
          "num_atoms": 25,
          "num_bonds": 24,
          "molecular_weight": 180.157
        }
      }
    },
    "processing_times": {
      "classical": { "117_ideal.mol2": 0.1234 },
      "quantum": { "117_ideal.mol2": 0.0567 }
    }
  }
}
```

## Frontend Integration

### JavaScript/React Example

```javascript
const uploadFiles = async (files, approach) => {
  const formData = new FormData();
  files.forEach((file) => formData.append("files", file));
  formData.append("approach", approach);

  try {
    const response = await fetch("http://localhost:8000/upload-and-process", {
      method: "POST",
      body: formData,
    });

    const result = await response.json();
    return result;
  } catch (error) {
    console.error("Upload error:", error);
  }
};

// Usage
const files = document.getElementById("fileInput").files;
const result = await uploadFiles(files, "both");
console.log("Session ID:", result.session_id);
```

## Testing

Run the test client to verify everything is working:

```cmd
python test_client.py
```

This will:

1. Test the health check endpoint
2. Upload the sample MOL2 file
3. Retrieve processing results
4. Generate and save summary plots

## File Structure

```
Backend/
├── main.py                 # FastAPI application
├── classical_processor.py  # Classical computing logic
├── quantum_processor.py    # Quantum computing logic
├── visualization_utils.py  # Plotting and visualization utilities
├── requirements.txt        # Python dependencies
├── start_server.bat       # Server startup script
├── test_client.py         # Test client for API
├── 117_ideal.mol2         # Sample MOL2 file
└── README.md              # This file
```

## Dependencies

- **FastAPI**: Web framework for building APIs
- **uvicorn**: ASGI server for running FastAPI
- **RDKit**: Chemistry toolkit for molecular processing
- **PennyLane**: Quantum computing framework
- **py3Dmol**: 3D molecular visualization
- **matplotlib**: Plotting and visualization
- **numpy**: Numerical computing

## Error Handling

The API includes comprehensive error handling:

- Invalid file formats are rejected
- Processing errors are captured and returned
- Session cleanup prevents memory leaks
- Detailed error messages for debugging

## Production Considerations

For production deployment:

1. **Security**: Configure CORS properly and add authentication
2. **Database**: Replace in-memory storage with a proper database
3. **File Storage**: Use cloud storage for uploaded files
4. **Logging**: Add comprehensive logging
5. **Monitoring**: Add health checks and monitoring
6. **Scaling**: Consider using multiple workers and load balancing

## License

This project is for educational and research purposes.
