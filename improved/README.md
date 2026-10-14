# CarViG

Integrated Navigation for INS/GNSS/Vision for Vehicle Positioning

## Overview

CarViG is an enhanced navigation framework that extends INS/GNSS fusion with visual information for vehicle positioning. It is designed for challenging driving environments where GNSS can be weak, intermittent, or degraded and where vision adds useful motion and pose constraints.

The project is a research-oriented extension of the IGNAV navigation stack, with added support for visual-assisted positioning and multi-sensor fusion.

## Features

- INS/GNSS loose-coupling and tight-coupling algorithms
- visual information-assisted positioning
- relative pose estimation for motion constraints
- OpenCV and Eigen integration
- robust filtering and smoothing for trajectory estimation
- modular design for experimentation and extension

## Requirements

- Linux environment
- CMake
- C++11
- OpenCV (2.4.3 or 3.x)
- Eigen3
- Python development libraries
- BLAS/LAPACK
- zlib, pthread, libpng, and JPEG support

## Quick Start

```bash
cd improved
mkdir -p build
cd build
cmake ..
make -j$(nproc)
```

The build produces executables in the `improved/bin/` directory.

## Executables

- `carvignav` — navigation executable for visual and inertial/GNSS fusion workflows
- `tc-postpos` — tightly-coupled post-processing executable

## Project Layout

```text
improved/
├── CMakeLists.txt
├── README.md
├── bin/
├── cmake/
├── data/
├── include/
├── src/
├── lib/
└── ...
```

## Typical Workflow

1. Prepare sensor data and configuration files.
2. Build the project.
3. Run the selected executable.
4. Analyze the resulting trajectory and pose estimates.

## Notes

- This project depends on OpenCV and Eigen for visual and numerical processing.
- It is intended for research, algorithm validation, and prototype navigation tasks.

## References

- RTKLIB: http://www.rtklib.com
- OpenCV: https://opencv.org/
- Eigen: https://eigen.tuxfamily.org/
- Libviso: http://www.cvlibs.net/software/libviso/

## Contact

- QQ: 1971129844
- Email: sujinglan@whu.edu.cn
