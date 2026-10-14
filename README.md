# IGNAV

Integrated Navigation for INS and GNSS

## Overview

IGNAV is a C/C++ library for integrated navigation and positioning based on INS (Inertial Navigation System) and GNSS measurements. It is designed for vehicle navigation and multi-sensor fusion applications, with emphasis on robust state estimation, filtering, and smoothing for high-precision positioning.

The project extends the traditional RTKLIB-based framework with additional navigation functions for inertial/GNSS integration, motion constraints, alignment initialization, and smoothing strategies. It provides a practical and extensible research-oriented software platform for INS/GNSS navigation algorithm development.

## Motivation

Accurate and continuous navigation is essential for intelligent vehicles, mobile robotics, and autonomous systems. In real-world scenarios, standalone GNSS may degrade in urban canyons, tunnels, and multipath-dominated environments, while inertial sensors can provide short-term continuity but drift over time. IGNAV addresses this challenge by combining inertial measurements and GNSS observations in a unified estimation framework.

## Features

### GNSS Processing

- Full support for multi-GNSS observations, including Beidou, GPS, GLONASS, Galileo, and QZSS
- SPP, PPP, DGPS, and RTK-based positioning modes
- Robust least-squares solutions for GNSS estimation
- Sequential outlier detection and repair strategies
- Clock bias and drift jump detection and compensation
- Adaptive stochastic error modeling for satellite measurements
- Adaptive Kalman filter noise estimation
- Motion-state discrimination, including static and dynamic conditions
- Reliable ambiguity search and fixing strategies
- Satellite observation quality monitoring
- Autonomous integrity monitoring for positioning quality assessment
- Parallel multi-stage fusion feedback positioning framework
- Forward/backward RTS smoothing for improved trajectory consistency

### INS/GNSS Integration

- INS/GNSS loose-coupling algorithms
- INS/GNSS tight-coupling algorithms, including SPP, PPP, DGPS, and RTK variants
- Robust Kalman filtering for GNSS-aided inertial navigation
- Motion constraints such as NHC, ZUPT, and ZARU
- Odometer-assisted navigation
- Magnetometer-assisted navigation
- Dual-antenna heading assistance
- Forward and backward INS mechanization
- Dynamic and static alignment initialization
- Forward/backward integrated filtering for INS/GNSS fusion
- RTS forward/backward smoothing
- Real-time integrated state estimation for position, velocity, and attitude

### Additional Capabilities

- Preliminary support for visual-information-assisted positioning and attitude determination
- Dynamic display of vehicle trajectories
- Modular design suitable for algorithm research and engineering prototyping

## Repository Layout

```text
.
├── CMakeLists.txt        # Build configuration
├── README.md             # Project documentation
├── include/              # Public headers and APIs
├── src/                  # Core source code
├── example/              # Example configurations, data, and demos
│   ├── conf/
│   ├── data/
│   └── demo/
├── bin/                  # Compiled executables and generated outputs
├── data/                 # Supporting data files
├── plot/                 # Plotting-related utilities
├── man/                  # Manual or documentation assets
├── improved/             # Enhanced or experimental implementations
├── lib/                  # Built library outputs
└── .gitignore
```

## Requirements

- Linux environment
- CMake
- BLAS, LAPACK, zlib, and pthread support
- Qt (required for plotting functionality)

## Build

From the repository root, configure and build the project with CMake:

```bash
mkdir -p build
cd build
cmake ..
make -j$(nproc)
```

The project builds the core library and executables, including:

- `navapp`
- `lc-rts`
- `lc-fbsm`

Executables are generated in the `bin/` directory by default.

## Quick Start

Before running processing, configure the relevant parameter files according to the RTKLIB format. Example configuration files are provided in the `example/conf` directory.

A typical execution example is:

```bash
./bin/navapp -o ../example/conf/navlib.conf -m 52716
```

The program uses the configuration file to initialize the navigation pipeline and processes the specified dataset. By default, outputs are stored in the same directory as the executable in `bin/`.

## Example Data and Demos

The project includes example data and demonstration resources in:

- `example/conf/`
- `example/data/`
- `example/demo/`

These resources provide a practical way to understand the data flow, configuration format, and expected output behavior.

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
