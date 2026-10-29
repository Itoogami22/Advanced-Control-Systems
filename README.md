# Advanced Control Systems

Coursework repository for the analysis of a 3-DOF PRR (Prismatic-Revolute-Revolute) robot manipulator. It contains MATLAB scripts, the robot's URDF model, figures, and a LaTeX report for each assignment.

## Repository contents

```text
Images/                       Figures used in the reports
MatlabScripts/                MATLAB implementations and robot URDF
Report/
  1_FirstAssignment/          Kinematics and Jacobians report
  2_SecondAssignment/         Energy and inertia report
  3_ThirdAssignment/          Equations of motion report
  4_FourthAssignment/         Dynamics and identification report
```

Each assignment has a matching `.m` script in `MatlabScripts/` and `.tex` report in its corresponding `Report/` folder:

1. **Kinematics:** Denavit-Hartenberg model, direct and inverse kinematics, and geometric and analytical Jacobians.
2. **Energy and inertia:** center-of-mass positions, partial Jacobians, kinetic and potential energy, and the generalized inertia matrix.
3. **Equations of motion:** dynamic matrices, equations of motion, and passivity/skew-symmetry analysis.
4. **Dynamics and identification:** Recursive Newton-Euler formulation, comparison with Euler-Lagrange dynamics, and least-squares parameter identification.

## Requirements

- MATLAB with **Robotics System Toolbox** and **Symbolic Math Toolbox**.
- A LaTeX distribution to compile the reports. The first report uses images from `Images/`.

## Run a MATLAB assignment

Set MATLAB's current folder to `MatlabScripts/`, then run the desired script:

```matlab
FirstAssignment
```

Use `SecondAssignment`, `ThirdAssignment`, or `FourthAssignment` for the other assignments. Run scripts from `MatlabScripts/` so MATLAB can locate the URDF model.

## Compile a report

Open the corresponding `.tex` file under `Report/` in your LaTeX editor and compile it with `pdflatex` (or an equivalent LaTeX engine).
