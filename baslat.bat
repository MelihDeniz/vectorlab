@echo off
cd /d "%~dp0"
title VectorLab
echo.
echo   VectorLab baslatiliyor...
echo   Tarayici acilmazsa: http://localhost:8501
echo   Kapatmak icin bu pencerede Ctrl+C
echo.
python -m streamlit run streamlit_app.py
if errorlevel 1 (
  echo.
  echo   Hata olustu. Bagimliliklar eksikse su komutu calistirin:
  echo   pip install -r requirements.txt
  echo.
  pause
)
