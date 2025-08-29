@echo off
title Qpharma - Master Startup Script
echo ================================================
echo          Qpharma - Quantum Drug Discovery
echo ================================================
echo.
echo This script will start both the backend and frontend servers.
echo Make sure you have Python 3.8+ and Node.js 16+ installed.
echo.
echo Backend will start at: http://localhost:8000
echo Frontend will start at: http://localhost:8080
echo.
echo Press any key to continue or Ctrl+C to cancel...
pause >nul

REM Start Backend Server in a new window
echo Starting Backend Server...
start "Qpharma Backend" cmd /c "cd Backend && start_server.bat"

REM Wait a few seconds for backend to start
echo Waiting for backend to initialize...
timeout /t 5 >nul

REM Start Frontend Server in a new window
echo Starting Frontend Server...
start "Qpharma Frontend" cmd /c "start_frontend.bat"

echo.
echo ================================================
echo Both servers are starting in separate windows...
echo ================================================
echo.
echo Backend: http://localhost:8000 (API Documentation: /docs)
echo Frontend: http://localhost:8080
echo.
echo To stop the servers, close both command prompt windows.
echo.
echo Happy coding! 🚀
echo.
pause
