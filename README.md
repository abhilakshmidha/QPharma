# Qpharma - Quantum Drug Discovery Platform

A full-stack application that combines quantum computing and classical methods for molecular analysis and drug discovery.

## Project Structure

```
Abhilaksh/
├── Backend/                 # FastAPI backend server
│   ├── main.py             # Main FastAPI application
│   ├── classical_processor.py
│   ├── quantum_processor.py
│   ├── visualization_utils.py
│   ├── requirements.txt
│   ├── start_server.bat
│   └── test_client.py
└── Frontend/               # React frontend application
    └── qubit-discovery-hub/
        ├── src/
        ├── package.json
        └── vite.config.ts
```

## Quick Start

### 1. Start the Backend Server

Navigate to the Backend directory and run:

```bash
cd Backend
start_server.bat
```

Or manually:

```bash
cd Backend
pip install -r requirements.txt
uvicorn main:app --host 0.0.0.0 --port 8000 --reload
```

The backend will be available at `http://localhost:8000`

### 2. Start the Frontend Development Server

Open a new terminal, navigate to the Frontend directory and run:

```bash
cd Frontend/qubit-discovery-hub
npm install
npm run dev
```

The frontend will be available at `http://localhost:8080`

## Features

### Backend (FastAPI)

- **File Upload**: Accept MOL2 molecular structure files
- **Classical Processing**: Traditional molecular simulation using RDKit
- **Quantum Processing**: Quantum simulation using PennyLane
- **3D Visualization**: Interactive molecular visualizations using py3Dmol
- **Performance Analysis**: Timing comparisons and statistical summaries
- **RESTful API**: JSON-based API with comprehensive error handling

### Frontend (React + TypeScript)

- **Modern UI**: Built with Tailwind CSS and shadcn/ui components
- **File Upload**: Drag & drop interface for MOL2 files
- **Analysis Options**: Choose between Classical, Quantum, or Both approaches
- **Real-time Results**: Live updates during processing
- **Interactive Visualizations**: 3D molecular structure viewers
- **Performance Charts**: Visual comparison of processing times
- **Responsive Design**: Works on desktop and mobile devices

## API Endpoints

### Core Endpoints

- `POST /upload-and-process` - Upload files and process with specified approach
- `GET /results/{session_id}` - Get detailed results for a processing session
- `GET /summary/{session_id}` - Get summary plots and statistics
- `DELETE /cleanup/{session_id}` - Clean up session data

### Specialized Endpoints

- `POST /process-classical` - Process files using classical approach only
- `POST /process-quantum` - Process files using quantum approach only
- `GET /health` - Health check endpoint

## Usage Workflow

1. **Upload MOL2 Files**: Drag and drop or select MOL2 molecular structure files
2. **Choose Approach**: Select Classical, Quantum, or Both processing methods
3. **View Results**: See processing times, molecular properties, and analysis results
4. **3D Visualization**: Open interactive 3D molecular structures in new windows
5. **Performance Analysis**: Generate summary charts comparing different approaches
6. **Download Results**: Export results as JSON files

## Technologies Used

### Backend

- **FastAPI**: Modern Python web framework
- **RDKit**: Chemistry toolkit for molecular processing
- **PennyLane**: Quantum machine learning library
- **py3Dmol**: 3D molecular visualization
- **Matplotlib**: Plotting and data visualization
- **Uvicorn**: ASGI server

### Frontend

- **React 18**: Modern React with hooks
- **TypeScript**: Type-safe JavaScript
- **Vite**: Fast build tool and dev server
- **Tailwind CSS**: Utility-first CSS framework
- **shadcn/ui**: Modern React component library
- **Lucide React**: Beautiful icons
- **React Query**: Server state management

## Development

### Backend Development

- Edit Python files in the `Backend/` directory
- The server runs with auto-reload enabled
- Test endpoints using the provided `test_client.py`
- API documentation available at `http://localhost:8000/docs`

### Frontend Development

- Edit React components in `Frontend/qubit-discovery-hub/src/`
- Hot module replacement (HMR) enabled
- TypeScript provides compile-time error checking
- Tailwind CSS for rapid UI development

### Testing

```bash
# Test backend API
cd Backend
python test_client.py

# Test frontend (run both servers first)
# Open http://localhost:8080 and upload a MOL2 file
```

## Configuration

### Environment Variables

No environment variables are required for basic setup. The application uses:

- Backend: `localhost:8000`
- Frontend: `localhost:8080`
- API proxy: `/api/*` routes to backend

### CORS Configuration

The backend is configured to allow requests from any origin for development.
For production, update the CORS configuration in `main.py`.

## Deployment

### Production Backend

```bash
# Install production ASGI server
pip install gunicorn

# Run with Gunicorn
gunicorn main:app -w 4 -k uvicorn.workers.UvicornWorker --bind 0.0.0.0:8000
```

### Production Frontend

```bash
# Build for production
npm run build

# Serve static files
# Deploy the `dist/` folder to your web server
```

## Troubleshooting

### Common Issues

1. **Port Already in Use**

   - Change ports in `start_server.bat` or `vite.config.ts`
   - Kill existing processes using the ports

2. **Python Package Installation Issues**

   - Ensure you have Python 3.8+ installed
   - Consider using a virtual environment
   - Install packages one by one if bulk installation fails

3. **Frontend Build Issues**

   - Ensure Node.js 16+ is installed
   - Clear npm cache: `npm cache clean --force`
   - Delete `node_modules` and run `npm install` again

4. **File Upload Issues**
   - Ensure MOL2 files are valid
   - Check file size limits (10MB default)
   - Verify backend server is running

### Support

- Check the console logs for detailed error messages
- Verify both servers are running before testing
- Use the test client to verify backend functionality

## License

This project is for educational and research purposes.
