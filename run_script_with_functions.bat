@echo off
echo ========================================
echo Starting Abaqus script in NoGUI mode...
echo ========================================
abaqus cae nogui=extract_example_with_functions.py > output.log 2>&1
echo ========================================
echo Script execution completed.
echo Log messages are stored in output.log.
pause