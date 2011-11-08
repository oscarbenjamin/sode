@echo off
REM Add to PATH and PYTHONPATH so that sode works when not installed
set PATH=%cd%\scripts;%PATH%
if "%PYTHONPATH%"=="" (
    set PYTHONPATH=%cd%
) else (
    set PYTHONPATH=%cd%;%PYTHONPATH%
)
