@echo off
echo Starting Qpharma Frontend Development Server...
echo ================================================

REM Navigate to the frontend directory
cd "Frontend\qubit-discovery-hub"

REM Check if Node.js is installed
node --version >nul 2>&1
if errorlevel 1 (
    echo Error: Node.js is not installed or not in PATH
    echo Please install Node.js 16+ from https://nodejs.org/
    pause
    exit /b 1
)

REM Check if npm is available
npm --version >nul 2>&1
if errorlevel 1 (
    echo Error: npm is not available
    pause
    exit /b 1
)

REM Install dependencies if node_modules doesn't exist
if not exist "node_modules" (
    echo Installing dependencies...
    npm install
)

REM Start the development server
echo Starting frontend development server...
echo Frontend will be available at: http://localhost:8080
echo Make sure the backend server is also running at: http://localhost:8000
echo Press Ctrl+C to stop the server
echo.

npm run dev

pause
