# Script to build necessary libraries and compile Backus-Gilbert scripts

# Generate shared object files for ZKCM interface
bash zcompile.sh

# Compile scripts
gcc "bgv6.c" -Wall -L. -linterface -lzcall -lstdc++ -lzkcm -lmpfr -lgmp -lgmpxx -o backus -fopenmp -lm -std=gnu99
gcc chread.c -o chread -lm
gcc covgen.c -o covgen -lm
gcc covgend.c -o covgend -lm
echo "done."
