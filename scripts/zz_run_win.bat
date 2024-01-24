echo on
copy ..\..\_build\tcsam02.exe
tcsam02 -rs -nox -configFile ..\MCI.inp -phase 5 -nohess -pin ..\tcsam02.pin
pause
