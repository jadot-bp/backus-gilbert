
## Authors

- [@jadot-bp](https://www.github.com/jadot-bp)
# backus-gilbert [hlt]

This repository is for the Backus-Gilbert code as used by FASTSUM for spectral reconstruction.


## Installation (Linux)

Clone the project

```bash
  git clone https://github.com/jadot-bp/backus-gilbert
```

Go to the project directory

```bash
  cd backus-gilbert
```

Switch to the `hlt` branch

```bash
  git checkout hlt
```

Download and install the following dependencies: 

*  GMP: https://gmplib.org/
* MPFR: https://www.mpfr.org/ 
* ZKCM: https://sourceforge.net/projects/zkcm/

Compile the base scripts along with the ZKCM connector

```bash
  bash compile.sh
  bash zcompile.sh
```
Done!

## Generating the kernel weight matrix:

This code does not perform the full HLT method, but instead performs a high-precision calculation of the inverse of the kernel weight matrix, which is offloaded to the file `out.kinv`.

At the minimum, the operation of the script requires an input file of the format

```bash
  Nt
  Ns
```

where `Nt` is the lattice temperature/lattice extent and `Ns` is the number of sample slices (The number of coefficients output will therefore be `Ns+1`).
If the Tikhonov whitening flag `tikh` is not set (equal to 1), the input file will also require the covariance matrix to be specified. The code assumes that the covariance matrix is supplied in the format:

First `Nt` lines are the main diagonal of the covariance matrix. The subsequent `(Nt+1)Nt/2` lines are the off-diagonals in row-major order (left-to-right).

## Common Issues:

If you encounter the error

```bash
  ./backus_lsq: error while loading shared libraries: libinterface.so: cannot open shared object file: No such file or directory
```

then you have not correctly set the library path. Check to see if `libinterface.so` is in the linked library path

```bash
  echo $LD_LIBRARY_PATH
```

If libinterface.so is not in the listed path or the listed path is empty, link the path using the following:

```bash
  cd backus-gilbert
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$PWD"
```

## Support

For support, email Ben Page at 900727@swansea.ac.uk.


