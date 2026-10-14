# IGNAV
Integrated Navigation for INS and GNSS

## Overview
IGNAV is an INS/GNSS integrated navigation algorithm library developed based on RTKLIB and implemented in C/C++. It is designed for vehicle navigation and positioning applications and supports a wide range of navigation filtering, initialization, and smoothing techniques.

The library focuses on integrating inertial sensors with GNSS measurements to improve navigation accuracy, robustness, and continuity in dynamic environments.

## Key Features
IGNAV includes the following major capabilities:

1. Loose-coupling algorithms
2. SPP, PPP, DGPS, and tightly coupled RTK algorithms
3. Odometer-assisted navigation
4. Magnetometer-assisted navigation
5. Motion constraints such as NHC, ZUPT, and ZARU
6. Doppler-aided INS/GNSS
7. Dual-antenna heading assistance
8. Static and dynamic alignment initialization
9. Forward and backward INS mechanization
10. Forward and backward INS/GNSS integrated filtering
11. Preliminary support for visual-information-assisted positioning and attitude determination
12. RTS forward/backward smoothing
13. Dynamic display of vehicle trajectory

## Repository Structure

- `src/` — core INS/GNSS algorithms and application sources
- `include/` — public headers and library interfaces
- `example/` — example configurations, datasets, and demo assets
- `bin/` — built executables and generated outputs
- `data/` — supporting data files
- `plot/` — plotting-related components
- `man/` — manual/documentation resources
- `improved/` — experimental or enhanced implementations
- `CMakeLists.txt` — build configuration for the project

## Requirements

- Linux environment
- CMake
- BLAS, LAPACK, zlib, and pthread support
- Qt (required for plot functionality)

The project follows the standard CMake build flow and is intended for Linux-based development and deployment.

## Build

From the repository root:

```bash
mkdir build
cd build
cmake ..
make -j$(nproc)
```

This will build the project libraries and executable tools, including:

- `navapp`
- `lc-rts`
- `lc-fbsm`

The binaries are output to the `bin/` directory by default.

## Configuration and Processing

Before running processing, configure the relevant parameter files according to the RTKLIB format. Example configuration files can be found in the `example/conf` directory.

By default, output results are written to the same directory as the executable in `bin/`.

## Run

Example command:

```bash
./bin/navapp -o ../example/conf/navlib.conf -m 52716
```

There are also demonstration videos and example assets under `example/demo`.

## References

[1] RTKLIB: http://www.rtklib.com

[2] P. D. Groves, Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems

[3] Wuhan MAP Navigation Technology Co., Ltd., http://www.whmpst.com/cn/

## Contact

Address: Satellite Navigation and Positioning Technology Research Center, Wuhan University

Researcher: Su Jinglan

QQ: 1971129844

Email: sujinglan@whu.edu.cn

## Note

If you need the professional version of IGNAV, please contact the author.
