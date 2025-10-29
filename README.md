# ckalman

The Kalman filter here is implemented in C style C++ for easy integration into C and C++ projects. 
Memory allocation is managed through `memory_t` and should give the user full control over memory usage.

Features:

- Self contained, including the math
  - LU decomposition based matrix inversion
- No dynamic memory allocation (user provides memory)

## State

- Alpha version; next step is to add unittests

