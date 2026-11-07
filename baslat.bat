@echo off
rem VectorLab'i kendi penceresinde acar (tarayici acilmaz).
rem Konsol gorunmesin diye pythonw kullanilir; sorun ciktiginda
rem hata kutusu belirir ve ayrinti data\desktop.log icine yazilir.
cd /d "%~dp0"
start "" pythonw.exe "%~dp0desktop.py"
