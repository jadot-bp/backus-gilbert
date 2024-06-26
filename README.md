
## Authors

- [@jadot-bp](https://www.github.com/jadot-bp)
# backus-gilbert

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

Download and install the following dependencies: 

*  GMP: https://gmplib.org/
* MPFR: https://www.mpfr.org/ 
* ZKCM: https://sourceforge.net/projects/zkcm/

Compile the base scripts along with the ZKCM connector using the Makefile

```bash
  make
```
Done!

## Common Issues:

If you encounter the error

```bash
  ./backus: error while loading shared libraries: libinterface.so: cannot open shared object file: No such file or directory
```

then you have not correctly set the library path. Check to see if `libinterface.so` is in the linked library path

```bash
  echo $LD_LIBRARY_PATH
```

If libinterface.so is not in the listed path or the listed path is empty, link the path using the following:

```bash
  cd backus-gilbert
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:lib/"
```

## Support

For support, email Ben Page at benjaminpage.acer@gmail.com.


