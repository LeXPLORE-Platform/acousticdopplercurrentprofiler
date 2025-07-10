:: UTF-8 characters (accents):
chcp 65001

@echo off
setlocal enabledelayedexpansion

:: Ensure correct location
:: cd "C:\Users\Seatronic 1147\Documents\Data_Lexplore\git\acoustic-doppler-current-profiler

:: Load input variables
call "scripts\input_batch_TD.bat"

:: Backup files
:: robocopy %in300% %raw300% *.LTA /E /MIR
:: robocopy %in600% %raw600% *.LTA /E /MIR

:: Process ADCP data
:: To process the most recent data only:
::%pythonenv% %script% live

:: To process everything:
%pythonenv% %script%
echo Data has been processed
pause

:: %pythonenv% %upload% -w

:: curl "https://api.datalakes-eawag.ch/update/599"
:: timeout 20
:: curl "https://api.datalakes-eawag.ch/update/600"
