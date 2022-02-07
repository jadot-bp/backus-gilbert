# Script to build necessary libraries and compile Backus-Gilbert scripts

# Generate shared object files for ZKCM interface
bash zcompile.sh

# Compile scripts
gcc "bgv6_leastsq.c" -Wall -L. -linterface -lzcall -lstdc++ -lzkcm -lmpfr -lgmp -lgmpxx -o backus_lsq -fopenmp -lm -std=gnu99
gcc "bgv6_spread.c" -Wall -L. -linterface -lzcall -lstdc++ -lzkcm -lmpfr -lgmp -lgmpxx -o backus_s -fopenmp -lm -std=gnu99
gcc chread.c -o chread -lm
gcc covgen.c -o covgen -lm
gcc covgend.c -o covgend -lm
echo "done."
