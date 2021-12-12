# Script to build necessary libraries and compile Backus-Gilbert scripts

# Generate shared object files for ZKCM interface
bash zcompile.sh


gcc "bgv5.c" -Wall -L. -linterface -lzcall -lstdc++ -lzkcm -lmpfr -lgmp -lgmpxx -o backus_v5 -fopenmp -lm -std=gnu99
gcc consol.c -o consol -lm
gcc bootstrap.c -o boots -lm
gcc bootstrapv2.c -o boots_v2 -lm
gcc compress.c -o compress -lm
echo "done."
