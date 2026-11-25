# Skid-Steer Simulator (MATLAB)

Interactive simulators for **four-wheel skid-steer robots**: compare kinematic
models, and explore a full-physics model where the instantaneous center of
rotation (ICR), wheel spin, skid, loads and motor currents *emerge* from friction,
terrain and drivetrain limits.

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=YOUR-USERNAME/skid-steer-sim)



## Try it in your browser (no install)

**[Open the web version](https://YOUR-USERNAME.github.io/skid-steer-sim/)**, or download
[`index.html`](index.html) and double-click it. It runs offline in any browser and
has the same physics as `skidSteerPro.m` (one view instead of two).

## Quick start (MATLAB)

```matlab
% open MATLAB in the folder containing these files, then run:
skidSteerPro        % full-physics interactive simulator (start here)
skidSteerSim        % kinematic model comparison (ideal / chi / slip / ICR)
skidSteerVisual     % ICR explorer: pick an ICR, see wheel rotations
```

**Requirements:** MATLAB R2017b or newer. No toolboxes. (Octave is not supported:
the GUIs use nested functions, `uitable` and `timer`.) No MATLAB installed? Use
**MATLAB Online** with the button above, or the Python reference `skid_steer_reference.py`.

## What each file does

| File | Purpose |
|---|---|
| `skidSteerPro.m` | Interactive time-domain simulator, full physics (below). |
| `skidSteerSim.m` | Scripted comparison of kinematic models on the same wheel commands, with plots, animation and a `chi` sweep. Edit `defaultParams()`. |
| `skidSteerVisual.m` | You choose the ICR (drag it); shows each wheel's rotation, skid and friction. Quasi-static, **educational**: it does not check the ICR is physically sustainable. |
| `index.html` | Single-file browser version of the full-physics simulator (also served by GitHub Pages). |
| `skid_steer_reference.py` | Headless NumPy port of the `skidSteerPro` physics for batch runs and cross-checks. |

## Physics in `skidSteerPro`

- Planar rigid body (vx, vy, yaw) with centripetal terms, CG offset, yaw-inertia factor, aero drag
- Dynamic wheel loads: longitudinal + lateral **load transfer** (CG height, front/rear split, suspension lag), wheel lift-off and static tip limits
- Slope angle and uphill direction
- **Anisotropic friction** ellipse (mu along / across the wheel) with a soil-shear model: friction builds with slip distance, giving real static friction
- **Terrain presets** (hard, gravel, sand, loam, snow): rolling resistance + Bekker sinkage/compaction resistance
- Contact-patch scrub moment vs wheel width
- **Drivetrain**: PI wheel-speed loop, DC motor (R, kt, back-EMF), gearbox (ratio, efficiency, reflected inertia), current limit, battery internal resistance
- Outputs: per-wheel rotation, current, torque, slip, load and forces; vehicle ICR; effective track factor chi; track ICRs; power breakdown

Frames: x forward, y left, yaw positive CCW. Wheel order FL, FR, RL, RR; wheel spin positive = rolling forward.

### Key equations

Kinematic ICR model (Mandow et al., 2007), with track speeds `vL`, `vR`:

```
w  = (vR - vL) / (yL - yR)
vx = (yL*vR - yR*vL) / (yL - yR)
vy = -xv * w
```
Ideal differential drive is `yL = b/2, yR = -b/2, xv = 0`. The effective track
factor is `chi = (yL - yR) / b`.

Wheel contact slip velocity and body dynamics (CG frame):

```
sx = (vx - w*yi) - r*omega_i        sy = vy + w*xi
m(dvx - w*vy) = sum(Fx)             m(dvy + w*vx) = sum(Fy)
Iz*dw = sum(xi*Fy - yi*Fx) + scrub moment
Jeff*d(omega_i) = tau_wheel - r*Fx_i - rolling-resistance torque
```

## Status and limitations

- The physics was cross-checked against a Python port (standstill, slope hold,
  load transfer, terrain drag, current-limit stall). The MATLAB GUIs have **not
  yet been tested across MATLAB versions**: please open an issue if you hit one.
- Motion is planar; pitch/roll dynamics and suspension travel are not simulated
  (load transfer is quasi-static with a first-order lag).
- No motor inductance, thermal effects, backlash or varying terrain.
- Defaults are plausible estimates, **not validated against a specific robot**.
  Terrain Bekker values come from published tables; tune mu, K and rolling
  resistance to your own drive tests.

## Fitting to your robot

Measure mass, CG height and position, wheel radius and width, motor kt and
resistance, gear ratio, battery voltage. Then drive a few turns at known wheel
speeds, read the measured yaw rate, and tune friction (`mu_x`, `mu_y`) and soil
shear `K` until the simulated effective `chi` and yaw rate match.

## References

- A. Mandow et al., "Experimental kinematics for wheeled skid-steer mobile robots", IROS 2007.
- J. Y. Wong, *Theory of Ground Vehicles*.
- M. G. Bekker, *Introduction to Terrain-Vehicle Systems*, 1969.
- Z. Janosi and B. Hanamoto, "The analytical determination of drawbar pull as a function of slip for tracked vehicles in deformable soils", 1961.

## License

MIT, see [LICENSE](LICENSE).
