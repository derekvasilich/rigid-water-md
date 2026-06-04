# rigid-water-md

A high-fidelity rigid-body molecular dynamics engine for liquid water.

This project implements a 4-site water model with explicit rotational dynamics using quaternions, periodic boundary conditions, neighbor cell decomposition, and configurable OpenGL visualization.

## Architecture Overview

- **Molecular representation**
  - Each water molecule is a rigid 4-site body (TIP4P-style): 1 oxygen, 2 hydrogens, and an ancillary interaction site.
  - Molecules contain both translational state (`cm`, `veloc`, `accel`) and full rotational state (`q`, `qv`, `qa`, `omega`, `inert`).
  - Rigid-body orientation is stored as a normalized quaternion, eliminating gimbal lock and supporting smooth rotation.

- **Physics model**
  - Configurable via `water.h` for either water or argon-like operation.
  - Supports Coulombic and Lennard-Jones-style interactions for water.
  - Periodic boundary conditions are enforced in all three dimensions.
  - A cell-list spatial decomposition reduces the force evaluation domain to local neighbor cells.

- **Numerical integration**
  - Translation uses a predictor-corrector integration scheme with support for both Verlet and Gear-style updates.
  - Rotation uses quaternion-based integration and normalization after each step.
  - Angular acceleration is computed from torque and then converted into body-fixed angular velocity.

- **Data and analytics**
  - Outputs statistical traces for energy, velocity, angular momentum, pressure/temperature, and radial distribution functions.
  - Supports saving and loading binary restart snapshots.
  - Built-in OpenGL rendering provides a real-time molecular view when enabled.

## Key Design Highlights

- **Quaternion-based rigid rotation**: Each molecule updates orientation via quaternion kinematics rather than Euler angles.
- **Rigid 4-site water model**: Simulates molecular geometry and internal site forces, not just point particles.
- **Neighbor grid / cell lists**: Uses a 3D cell array to organize molecules for efficient pair interactions.
- **Periodic domain wrapping**: Molecule positions are wrapped across boundaries automatically.
- **Flexible runtime mode**: Command-line selectable display rate, simulation length, statistics mode, and optional OpenGL output.

## Build Instructions

### macOS / Linux

Required libraries:
- `GLUT`
- `OpenGL`

Build with:

```sh
make
```

This produces the executable `water.exe`.

### Windows

A `make.bat` stub is included for Windows-style builds.

## Usage

```sh
./water <dispStep> <totalStep> [statsType] [useGL] [restartFile]
```

Example:

```sh
./water 1000 100000 1 1
```

Parameters:
- `dispStep`: number of timesteps between redraws or progress updates
- `totalStep`: total number of simulation timesteps
- `statsType`: optional bitmask for statistics output
- `useGL`: `0` or `1` to disable/enable OpenGL rendering
- `restartFile`: optional binary snapshot file to resume a simulation

## Runtime Output

The simulation generates statistics logs for:
- total energy
- center-of-mass velocities
- momentum
- angular velocity
- pressure and temperature
- radial distribution functions (RDF)

## Configuration Points

`water.h` controls the simulation behavior:
- `NUM_MOLECULES`: system size
- `DENSITY`, `TEMPERATURE`: thermodynamic conditions
- `ROTATIONAL_MOTION`: enable/disable rotational dynamics
- `COULOMB_POTENTIALS`: enable electrostatic interactions
- `CUTOFF_DISTANCE`: interaction cutoff radius
- `DT`: integration timestep

## Why this project is architecturally strong

- Reflection of a mature MD architecture with separate concerns for geometry, forces, integration, and output.
- A carefully structured rigid body representation that tracks both center-of-mass motion and rotational kinetics.
- A performance-aware force loop using spatial cells and periodic neighbor scanning.
- A small but complete research-style prototype with restart capability and configurable analytics.

## Notes

- The code is intentionally low-level C and optimized around explicit structures and arrays.
- `water.c` also includes a legacy OpenGL front-end, turning the simulation into a live visualization tool.
- The project is a great foundation for extending to faster neighbor lists, Ewald summation, or GPU-accelerated force computation.

## History and Future Improvements

This project was designed as an early research-quality MD prototype that blends physics fidelity with a compact C implementation.

- **Historical context**: written when quaternion-based rigid-body integration was less common in lightweight MD demos, making it a strong example of realistic molecular simulation for its time.
- **Evolution potential**:
  - add a modern Verlet/Leapfrog integrator for rotation using quaternions with constraint stabilization
  - replace the brute-force cell list with linked-cell or Verlet neighbor lists for larger system scalability
  - integrate long-range electrostatics (Ewald, PPPM, or PME) for accurate water properties
  - migrate rendering to a modern OpenGL/GLFW or WebGL frontend for interactive visualization
  - add unit tests and configure an automated build pipeline for reproducibility

This makes the project both a capable legacy MD engine and an excellent starting point for a modern simulation research rewrite.
