"%PYTHON%" setup.py install
if errorlevel 1 exit 1

ln -s $PREFIX/CG2AT $PREFIX/test_sym
if errorlevel 1 exit 1

echo $PREFIX/AAAAAAAAAAAAAAAAAAAAAAAAA
if errorlevel 1 exit 1
