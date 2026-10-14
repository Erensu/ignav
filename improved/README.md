# CarViG

Integrated Navigation for INS/GNSS/Vision for Vehicle Positioning

## Overview

CarViG is an enhanced navigation and positioning framework that extends inertial navigation and GNSS fusion with visual information for vehicle localization. It is designed for high-precision positioning in complex driving environments where GNSS signals may be weak, intermittent, or degraded, and where vision provides complementary constraints for motion and pose estimation.

The project integrates:

- INS/GNSS fusion for robust navigation
- vision-based relative pose estimation
- multi-sensor state estimation for vehicle positioning
- smoothing and filtering strategies for trajectory consistency

CarViG is a research-oriented extension of the IGNAV navigation framework, with additional support for visual-assisted positioning and advanced perception-driven fusion.

## Motivation

Vehicle positioning systems must remain accurate and continuous in challenging conditions such as urban canyons, tunnels, dense foliage, and dynamic motion. GNSS alone is vulnerable to multipath and signal obstruction, while inertial sensors suffer from drift over time. By incorporating visual observations and relative pose constraints, CarViG improves continuity, robustness, and positioning accuracy for autonomous and intelligent vehicle applications.

## Key Features

### INS/GNSS Fusion

- loose-coupling and tight-coupling INS/GNSS algorithms
- GNSS-based positioning support within the navigation pipeline
- robust filtering for sensor fusion and anomaly handling
- motion constraint support, including NHC, ZUPT, and ZARU related strategies
- forward and backward mechanization workflows
- dynamic alignment initialization and fusion-based smoothing

### Visual-Aided Navigation

- visual information-assisted positioning
- relative pose estimation for motion constraints
- integration with OpenCV-based computer vision pipelines
- support for feature-based motion estimation and geometric constraint recovery
- compatibility with visual odometry / stereo-related positioning workflows

### Scientific and Engineering Components

- CMake-based build system
- OpenCV and Eigen integration
- support for Python-enabled components and numerical libraries
- modular library design for research extension
- prebuilt support for custom relative-pose and visual library modules

## Architecture

The project structure includes several major components:

```text
improved/
├── CMakeLists.txt
├── README.md
├── bin/                  # built executables
├── cmake/                # CMake helper modules
├── data/                 # example and support data
├── include/              # project headers
├── src/
│   ├── ins-gnss/         # INS/GNSS algorithms and application code
│   ├── libviso/          # libviso integration
│   ├── relative-pose/    # relative pose estimation module
│   └── ...
├── lib/                  # output libraries
└── ...
```

## Dependencies

CarViG relies on the following core libraries and toolchains:

- CMake
- C++11
- OpenCV (version 2.4.3 or 3.x)
- Eigen3
- Python development libraries
- BLAS/LAPACK
- zlib, pthread, and related system libraries
- libpng / JPEG support for visual processing

## Build

From the project root (`improved/`):

```bash
mkdir -p build
cd build
cmake ..
make -j$(nproc)
```

The build system compiles the main library and associated executables, including:

- `carvignav`
- `tc-postpos`

The output binaries are stored in the `bin/` directory by default.

## Main Executables

### `carvignav`

A navigation executable designed for visual-inertial/GNSS integrated positioning workflows.

### `tc-postpos`

A tightly-coupled post-processing executable for improved positioning results based on multi-sensor integration.

## Example Workflow

The project is intended for research and prototype evaluation in vehicle navigation tasks. Typical use involves:

1. preparing configuration files and sensor data
2. constructing the navigation system state
3. running the visual and inertial/GNSS fusion pipeline
4. analyzing the generated trajectory and pose estimates

## References

This project builds upon established tools and methods in navigation and computer vision, including:

- RTKLIB: http://www.rtklib.com
- OpenCV: https://opencv.org/
- Eigen: https://eigen.tuxfamily.org/
- Libviso: http://www.cvlibs.net/software/libviso/

## Contact

The project is associated with the navigation research work of Su Jinglan.

QQ: 1971129844

Email: sujinglan@whu.edu.cn

## Note

CarViG is an extended research-oriented version of the IGNAV navigation framework, emphasizing visual enhancement and multi-sensor fusion for vehicle positioning.
