# Ai Thinker RD03D mmWave Radar Module

A Kalman filter optimizes raw data from your RD03D mmWave sensor to provide smooth, highly accurate, and lag-free tracking.

## Why You Need It

- **Removes Noise**: mmWave sensors naturally suffer from environmental clutter, thermal noise, and signal fluctuations.
- **Prevents Jumps**: Stops the tracked target from "teleporting" or jittering between frames.
- **Handles Dropouts**: If a wall or object briefly blocks the target, the filter estimates the target's position until the signal returns.
- **Provides Velocity**: Calculates how fast an object is moving, even if the sensor only reports distance.

## Specific Use Cases

### 🎯 Micro-Movement and Presence Detection

The RD03D detects breathing and minor shifts. A Kalman filter separates these tiny, rhythmic human movements from random environmental static, preventing false negatives (for example, thinking a room is empty when someone is sleeping).

### 🚷 Zone and Boundary Tracking

If you divide a room into automation zones (for example, turning on lights only when near the desk), raw sensor drift can cause lights to flicker on and off. The filter smooths boundary lines so transitions are clean.

### 🏃‍♂️ Multi-Target Trajectory Smoothing

When tracking people moving through a space, raw coordinates can jump around. The filter creates a smooth, logical path (trajectory) for each person, making automation triggers highly reliable.

