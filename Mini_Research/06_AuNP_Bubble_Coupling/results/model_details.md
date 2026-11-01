# Thermal model, derivation and SI interface

## 1. Geometry and discretization

The bubble center is physically at zero. A spherical particle of radius a has
center radius r_np=R+a+d, where d is always the **surface gap**. In the numerical
meridian, the bubble center is (0,-R), its near apex (0,0), and the particle center
(0,a+d). Rotating the meridian around its cylindrical axis includes the full
particle and bubble. The liquid cylinder is bounded by radial coordinate R+pad,
bottom -2R-pad and top 2a+d+pad. Production pad=1 µm. A Dirichlet outer bath has
zero perturbation temperature; insulation is a separate domain control.

Triangle meshes polygonal sphere boundaries and seeded logarithmic rings, with
`pq20Q`. These are finite chordal spheres, not analytic curved elements. Increasing
the refinement factor changes all local arcs and radial seeds. Each triangle has
linear basis N_i, area A and cylindrical radial coordinates rho_i. The exact
diffusion matrix and row-sum lumped water capacity are

`K_ij = k*2*pi*A*mean(rho)*(grad(N_i) dot grad(N_j))`,
`C_i += rho_l*cp_l*2*pi*A*(sum(rho_j)+rho_i)/12`.

Boundary edge length l contributes nodal area `w_i=2*pi*l*(2*rho_i+rho_j)/6`.
Au links are G*w_i between each boundary water node and a single gold node with
capacity 4*pi*a³*rho_Au*cp_Au/3. The same link enters both equations with opposite
sign. Bubble links h*w_i connect water to a fixed zero gas perturbation. Eliminating
fixed outer bath nodes supplies a separately calculated bath-flux vector. It is
not defined by energy conservation. Summing all active equations gives the
independent input/storage/gas/bath balance. No ALE work or moving-water energy
is implied: the domain is fixed throughout each response.

P1 FEM has tiny nonphysical signed far-tail undershoots. All are saved, never
clipped. Refined meshes reduce them. Early exponentially small far-source tails
are not relatively resolved; they are not interpreted as measurable energy.

## 2. Heat source and signed gas closure

Read-only MR5 imports supply actual `optical_cross_sections(radius_m,
wavelength_m,n_medium=1.33,model='mie')`, `absorbed_power`, `pulse_energy`,
`ThermalProperties`, `thermal_scales` and the isolated `RadialHeatModel`.
Each particle's local intensity is specified after applicable losses. No
macroscopic layer energy is used. `E_abs=F_local*sigma_abs`; the reference
cross section is 3.4337324474063743e-15 m². The thermal top hat is applied exactly
on [0,duration], with a separate integration after pulse end.

The liquid calculation is **laser tagged**: all initial perturbations are zero,
and gas perturbation is zero at the Robin sink. Positive liquid-to-gas power is
`delta_q = integral(h*theta_water dA)` with n_g from gas to water. The absolute
physical law would be `h*(T_water-T_g)` and can be negative. The implemented gas
diagnostic uses the controlled small-increment option:

`q_total = -4*pi*k_g*R*(T_g-T0) + delta_q`.

The baseline is the inherited signed MR4 closure, **not** the absolute solution
of a new hot-gas liquid BVP. The full Robin flux is not added again. At nominal
h=k_g/R, a uniform stationary exterior liquid resistance changes the baseline
conductance by the factor 1/(1+k_g/k_l)=0.9578544061. The 4.2146% difference is
an analytic matching diagnostic, not a bound on the unsteady error. The unknown
gas boundary layer and hot transport make h a model sensitivity.

## 3. Particle motion and frozen-geometry admission

For passive incompressible flow, 3*r²*rdot=3*R²*Rdot, so r³-R³=constant and
`gap_dot=((R/r)**2-1)*Rdot`. Exterior gap therefore increases during collapse.
Fixed laboratory center and prescribed capture `r=R+a+d0` are distinct.
Capture has no derived force law here. All contacts/intersections stop the
liquid-only snapshot branch and receive their own status. No particle is moved
through gas while retaining an Au-water model.

Event snapshots freeze R and gap at the pulse center. They do **not** solve
the moving-boundary equation. Event-aware saved nodes are used for Pe and
radius/gap displacement, including nanosecond extrema. Gates: Pe<0.1,
relative radius change<0.01 and relative gap change<0.1 over 2 microseconds;
1.3 bar remains formal regardless. Relative radius change over diffusion time
is recorded separately; NaN means that time exceeds the saved domain. None
of the late strong-collapse traces satisfies a physical frozen-flow admission.

## 4. First-order gas diagnostic

The saturated fixed NASA9 air/H2O inventory is read from MR4. Total initial
pressure is p0+2*sigma/R0, and no additional constant vapor partial pressure
is added. A no-laser window starts at tc-3 µs using a new PCHIP interpolation
of the saved MR4 initial temperature, then integrates the common baseline on
the published MR2 Hermite mechanics. The window requests tc+2 µs. Its output
is not called a reproduction of the original dense polynomials or a five-cycle
global maximum.

With `f(T)=N(T)/C(T)`, `N=-3*n*Ru*T*U/R-K*(T-T0)` and K=4*pi*k_g*R,
`f_T=N'/C-N*C'/C²`. The exact NASA polynomial derivative supplies C'.
The first-order increment obeys `delta_T_dot=f_T(T_base)*delta_T+q_laser/C`.
This preserves the scalar comparison principle for positive sources and avoids
subtracting noisy large temperatures. Incremental work and loss use
`3*n*Ru*(U/R)*delta_T` and `K*delta_T`, with sensible increment `C(T_base)*delta_T`.
The independently accumulated ledger checks their first-order first law.
This does not include gas perturbation feedback into the liquid PDE.

A is a **transient early-heating diagnostic**, not a new equilibrium transformation
of R0, vapor pressure or gas composition. A/B use retained diffusion histories
at early/near/late pulse centers. C uses
a finite-interval normalized Gaussian with either the computed B energy over
5 µs or the ideal local absorbed energy. No energy is tuned to a temperature
target. The upper NASA boundary applies to T_base+delta_T; actual stopping
time and censoring temperature are saved. No missing tail is extrapolated.
Large pressure mismatch with the mechanical polytropic gas is retained.

## 5. API and memory

`FrozenHeatModel(Geometry(...), h_W_m2K, props, refinement, ...)` builds the
spatial state/operator. `pulse(E_abs_J, duration_s, time_s, ...)` evolves it
from zero perturbation; returns temperature rises K, powers W, signed energies J.
With `save_field=True`, it also returns all spatial nodes, triangles, boundary
markers, nodal areas and liquid temperatures. Boundary markers: bubble=1,
bath=2, axis=3, Au=4. Storage on omitted bath nodes is excluded consistently
from the active semidiscrete domain; domain refinement tests that approximation.

`FrozenHeatModel.from_saved_field(path, geometry=...)` reassembles the
stationary operator directly from archived nodes, triangles, edges and markers.
The original reference file requires explicit geometry; other fields carry it
in metadata. The saved mesh is a numerical input. A remesher is not called and
freshly generated node lists are not required to match between platforms. The
field test independently integrates water energy with degree-2 triangle
quadrature; three archived-mesh replays compare observable max norms at 2e-5
relative tolerance, rather than claiming remesher bitwise identity.

`ThermalHistory(time_s, laser_W, baseline_W, metadata,
reference_temperature_K=None, reference_radius_m=None)` stores thermal memory
on a closed interval. `baseline_W` is a **stored reference-state trace**.
The v2 return/NPZ fields are `q_wall_reference_baseline_W`,
`q_wall_reference_sum_W`, `q_wall_laser_increment_W` and reference-only
positive/negative parts. The reference sum is

`q_ref(t) = -4*pi*k_g*R_ref(t)*(T_ref(t)-T0) + delta_q_frozen(t)`.

It is not a physical total flux evaluated at a future supplied Tg. New exports
contain no `q_wall_total_W`. `current_state_total_flux_available` is always false.
In the retained original diagnostics, T_ref is the paired exported Tb+delta_T
and R_ref is that trajectory's prescribed R(t), both now explicitly archived;
they are distinct from the thermal solver's frozen geometry. No numerical source
trace is changed by the renaming. Older field names can be loaded but are always
returned as reference diagnostics with the v2 contract.

For an ODE evaluating its current-state baseline loss, use only

```python
qnet_W = -heat_loss_at_current_state_W + history.laser_increment_W(t_s)
```

The stored baseline is not added again. The changed-gas-temperature test checks
this rule. Neither the scalar increment nor the reference sum includes the
omitted gas-temperature response of the liquid. PCHIP, no extrapolation and
frozen-radius checks remain in force. Future moving-boundary coupling requires
MR7; this interface is not a universal dynamic gas/liquid source.

## 6. Saved-trace estimate of omitted liquid feedback

For a stationary unobstructed sphere in infinite water, uniform prescribed gas
perturbation has transfer function h/[h+k/R+k*sqrt(s/alpha)] from gas to surface
water. Its step response is

`S(t)=h/(h+k/R)*(1-erfcx((1/R+h/k)*sqrt(alpha*t)))`.

The surface rise is the Stieltjes convolution of S with saved delta_Tg. Piecewise
linear reconstruction and four-point Gauss quadrature evaluate it at two time
resolutions. This is a one-pass analytic diagnostic with prescribed input, no gas
update, no new finite-particle PDE, no parameter optimization and no MR7 solution.
It excludes finite-Au obstruction and finite-domain/motion effects.

The extra positive water-return term is H*theta_surface, H=4*pi*R²*h. The gas ODE already includes its own
negative -4*pi*k_g*R(t)*delta_Tg loss derivative; no second direct gas-loss term
is added in this estimate. `gas_feedback_estimate.csv` separately records
the already counted loss from 4*pi*k_g*R(t)*delta_Tg. For B, estimated integrated
return/source energy is about 0.2323%,0.09284%,0.000225% at 1.0,1.15,1.3 bar,
respectively. The last value covers only the censored observed interval.
A maximum-principle bound H*integral(max_so_far(delta_Tg)) is much looser: about
74.3%,12.1%,0.127%. Those bounds apply to the prescribed-trace stationary control;
the small analytic estimates are not rigorous error bounds for changing geometry,
and no universal negligibility claim follows. Original gas-temperature histories
remain unchanged.
