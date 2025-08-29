@echo off
echo Starting Molecular Simulation FastAPI Server...
echo ================================================

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo Error: Python is not installed or not in PATH
    pause
    exit /b 1
)

REM Install requirements if they don't exist
echo Installing/checking requirements...
pip install -r requirements.txt

REM Start the FastAPI server
echo Starting server on http://localhost:8000
echo Press Ctrl+C to stop the server
echo.

uvicorn main:app --host 0.0.0.0 --port 8000 --reload

pause
