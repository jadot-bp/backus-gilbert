# Script for creating shared object files for the ZKCM interface

# Compile ZKCM call shared object file
g++ zcall.cpp -lzkcm -lmpfr -lgmp -lgmpxx -o libzcall.so -fPIC -shared

# Compile C to C++ interface shared object file
g++ interface.cpp -o libinterface.so -L. -lzcall -fPIC -shared
