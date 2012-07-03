@echo off
REM copyright Oscar Benjamin 2011
REM
REM run sode as a python script.
REM
REM Because of:
REM STDIN/STDOUT Redirection May Not Work If Started from a File Association
REM http://support.microsoft.com/kb/321788
REM http://bugs.python.org/issue1012692
python -m sode.examples.pyexamples %*
