# Mathematical decisions recorded before implementation

The accepted starting snapshot is 493dbdc213f3789b7195d93760bd5c0c867cb98c;
local HEAD, origin/main and the independently read GitHub main agree. The full
tracked-file inventory is byte checked against this snapshot. No upstream writes
or parameter fitting are permitted. This is a reduced, conditional model.

## Unknowns and conservation

Level 3 integrates R,U,Tg, liquid cell temperatures and one representative Au
temperature for N identical, angularly separated particles. Local N is a count.
Level 2 uses the same no-laser Level 3 R,U history; both thermal directions remain.
The thermal energy balance includes prescribed mechanical work, not the entire
driven fluid's kinetic energy. No retaining-force work is modeled.

With u=R^2 U/r^2 and xi=(r^3-R^3)/3, Dxi/Dt=r^2 u-R^2 U=0.
Since dxi/dr=r^2, d/dr=r^2 d/dxi, and the material time derivative removes
advection, rho cp DtT=(k/r^2) d_r(r^2 d_r T) becomes
rho cp partial_t T|xi=partial_xi(k r^4 partial_xi T).
xi has units m^3; k r^4 T_xi has units W; its xi derivative is W/m^3.
Material volume is 4 pi dxi and cell heat capacity is 4 pi rho cp Delta xi.
All capacities are positive and independent of R. Every interior face transfers
the same power with opposite signs. Exact spherical thermal resistances between
cell centers are used; their current radii change at every RHS evaluation.

At the inner face, Qlg is positive into gas. The gas conductance H=4 pi eta kg R
is in series with the half-cell liquid resistance. Tw=Tg+Qlg/H is the reconstructed
surface temperature. The sole gas heat flux is Qlg: MR4 heat_loss and MR6 stored
reference baselines/sums are never added. The outer fixed material face has an
explicit conduction resistance to T0; Qbath is positive out of water. Gas pressure
is n Ru Tg/V, with all fixed H2O included and no separate pv. Log differentiation
gives pdot=pg(Tdot/Tg-3U/R), used in the full viscous Keller--Miksis equation.

## Particle reduction and motion

The finite Au surface is averaged over concentric shells around the bubble.
For a center rnp and radius a its radial surface-weight density is r/(2 rnp a)
on [rnp-a,rnp+a]. Cell weights integrate this density against piecewise linear
radial water reconstruction, constant between each boundary and nearest center,
and sum to one. The adjoint deposition uses the identical weights. This smooths
the derivative jumps of piecewise constant reconstruction at cell edges.
Tbar=sum(w_i Tl_i), QAu=4 pi a^2 G(TAu-Tbar); identical N multiplies both
the water source and total Au storage, once. Water occupies the entire shell;
the tiny particle's excluded volume and angular thermal resistance are omitted.
This can severely change local Au cooling and near-contact transfer. Compare
with archived MR6 fields at matched parameters without fitting weights or G.
The result does not resolve a local nanoparticle hotspot or its nucleation risk.

Passive rnp^3-R^3=constant, fixed laboratory center, and prescribed rnp=R+a+d0
are distinct. Fixed/capture scenarios have unverified external forces. Surface
gap <=0 or outer source escape terminates the branch; no gas-host continuation.
N>1 denotes a ring normal to the beam, with separation, projected shadowing and
single-pass extinction/beam checks. Isolated-sphere optics remain conditional.

## Reuse, hierarchy and screens

Read-only adapters reuse MR4 Mixture (NASA9 coefficients, energy, EOS and domains),
MR4 original solve/diagnostics for Level 0 regression, and MR5 optical_cross_sections,
absorbed_power, pulse_energy, ThermalProperties, thermal_scales. MR6 analytic point
controls, saved fields and source contract v2 are reference checks only. Level 1
is retained as the upstream diagnostic, not relabeled as an MR7 moving solution.

Wet gas stops at its 6000 K coefficient boundary; dry air stops at 20000 K;
ground-state dry Ar has constant Cv but no chemical validity at high temperature.
The lower coefficient boundary also stops. Accepted states never extrapolate.
Gold moderate rise, bulk melting reference, water temperature reference, isolated
optics, low liquid Mach, fixed composition, one-temperature pulse duration and
spatial averaging are separate screens. Material violations terminate the primary
branch; explicitly named formal continuations may cross reference screens with
fixed hypothetical properties, while contact/domain events always terminate.

Five cycles are requested. A censored pair has undefined global Delta Tmax;
only a shared-prefix diagnostic is allowed. Pulse center is fixed in absolute
time from the no-laser event and is never recentered on heated mechanics.
Weak-source results use a directional complex-step tangent about the no-laser
solution plus independent finite-fluence differences. This is a first-order
sensitivity, explicitly distinct from a finite nonlinear paired result.

No chemistry, CFD, MR8, broad optimization or experimental calibration is added.

## Strong-collapse resolution finding

The 1.15 bar nominal early-pulse derivative changes by about 0.6 percent from
48 to 96 material cells and retains its sign in an independent Radau calculation.
The extremely weak pulse-centered-at-collapse tail changes by nearly an order
of magnitude on the same refinement. That late tail is explicitly unresolved;
it is not promoted to a quantitative positive enhancement. The first water cell
expands in physical thickness during collapse despite its fixed material volume.
Short diffusion distances therefore require a separate collapse-time mesh test.

For censored gas histories, a variational temperature derivative at the fixed
reference endpoint is a formal endpoint diagnostic. It is not the derivative of
a five-cycle maximum, nor the derivative of the coefficient-censored maximum
(which is fixed at the boundary value). Global derivatives remain unavailable.
