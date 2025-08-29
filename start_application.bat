@echo off
echo.
echo ==============================================
echo  Quantum Drug Discovery Platform - Startup
echo ==============================================
echo.

echo Starting Backend Server (FastAPI)...
cd /d "%~dp0Backend"
start "Backend Server" cmd /k "python main.py"

echo Waiting for backend to initialize...
timeout /t 3 /nobreak >nul

echo.
echo Starting Frontend Server (React + Vite)...
cd /d "%~dp0Frontend\qubit-discovery-hub"
start "Frontend Server" cmd /k "npm run dev"

echo.
echo ==============================================
echo  Both servers are starting up...
echo  
echo  Backend:  http://localhost:8000
echo  Frontend: http://localhost:8080
echo  API Docs: http://localhost:8000/docs
echo.
echo  The application will be ready in ~10 seconds
echo ==============================================
echo.

timeout /t 5 /nobreak >nul

echo Opening application...
start "" "http://localhost:8080"

echo.
echo Setup complete! 
echo Press any key to exit this window...
pause >nul
