# IGNAV

Integrated Navigation for INS and GNSS

## Overview

IGNAV is a C/C++ library for INS/GNSS integrated navigation and positioning. It targets vehicle navigation and multi-sensor fusion applications, with emphasis on robust state estimation, filtering, and smoothing for high-precision positioning.

The project is built around RTKLIB-inspired GNSS processing and extends it with INS/GNSS fusion, motion constraints, alignment initialization, and smoothing methods.

## Features

- INS/GNSS loose-coupling and tight-coupling algorithms
- SPP, PPP, DGPS, and RTK-based processing
- Odometer-assisted navigation
- Magnetometer-assisted navigation
- Motion constraints including NHC, ZUPT, and ZARU
- Doppler-aided INS/GNSS fusion
- Dual-antenna heading assistance
- Static and dynamic alignment initialization
- Forward and backward INS mechanization
- Forward/backward integrated filtering and RTS smoothing
- Visual-information-assisted positioning support
- Dynamic trajectory display

## Requirements

- Linux environment
- CMake
- BLAS, LAPACK, zlib, and pthread support
- Qt for plotting functionality

## Quick Start

```bash
mkdir -p build
cd build
cmake ..
make -j$(nproc)
```

After building, the executables are generated in the `bin/` directory.

## Run Example

```bash
./bin/navapp -o ../example/conf/navlib.conf -m 52716
```

This command runs the navigation application using the example configuration under `example/conf`.

## Project Layout

```text
.
├── CMakeLists.txt
├── README.md
├── include/
├── src/
├── example/
├── bin/
├── data/
├── plot/
├── man/
├── improved/
└── lib/
```

## Example Data and Demos

- `example/conf/` — configuration files
- `example/data/` — sample datasets
- `example/demo/` — demo videos and example assets

## Notes

- Output files are written to the same directory as the executable in `bin/` by default.
- The project is intended for Linux-based development and deployment.

## References

- RTKLIB: http://www.rtklib.com
- P. D. Groves, Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems
- Wuhan MAP Navigation Technology Co., Ltd.: http://www.whmpst.com/cn/

## Contact

- Researcher: Su Jinglan
- QQ: 1971129844
- Email: sujinglan@whu.edu.cn
