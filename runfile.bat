@echo off
setlocal enabledelayedexpansion

:: Ensure correct location
cd "C:\Users\Seatronic 1147\Documents\Data_Lexplore\git\acoustic-doppler-current-profiler

:: Load input variables
call "C:\Users\Seatronic 1147\Documents\Data_Lexplore\git\acoustic-doppler-current-profiler\scripts\input_batch.bat"

:: Backup files
robocopy %in300% %raw300% *.LTA /E /MIR
robocopy %in600% %raw600% *.LTA /E /MIR

:: Process ADCP data
%pythonenv% %script% live

%pythonenv% %upload% -w

curl "https://api.datalakes-eawag.ch/update/599"
timeout 20
curl "https://api.datalakes-eawag.ch/update/600"
