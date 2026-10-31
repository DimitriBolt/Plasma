"""Generate data-derived numerical audit, trajectory admission and handoff facts."""
import json
import numpy as np
from .model import ROOT, Parameters
from .experiments import RESULTS, read_csv, write_csv, save_json


def main():
    rows=read_csv(RESULTS/'amplitude_summary.csv')
    cv=read_csv(RESULTS/'convergence.csv')
    admitted=[]
    anchors={float(r['Pa_bar']) for r in cv}
    for r in rows:
        complete=r['completed']=='True'; bar=float(r['Pa_bar']); model=r['model']
        if not complete:
            use='diagnostic_partial_interval_only'; reason='Requested five-cycle mechanics incomplete; no global minimum or full input'
        elif model=='RP':
            use='reference_equation_only'; reason='Incompressible reference; physical applicability must be assessed independently'
        elif bar<=1.:
            use='conditional_baseline'; reason='Small computed compressibility diagnostics; gas, spherical and transport assumptions remain'
        elif bar<=1.15:
            use='conditional_sensitivity'; reason='Transition scenario; event-specific closure and liquid-model comparison required'
        else:
            use='formal_extreme_scenario'; reason='Increasing compressibility diagnostics; Gilmore/EOS and gas-closure checks required before physical interpretation'
        admitted.append(dict(model=model,Pa_bar=bar,completed=complete,
            numerical_status='independent_integrator_anchor' if bar in anchors or (model=='RP' and bar==1.295) else 'production_run_only',
            max_wall_Mach=r['max_wall_Mach'],pressure_bulk_modulus_proxy=r['max_liquid_pressure_change_over_bulk_modulus'],
            equation_output='Observed t,R,U and refined extrema retained; global only if completed',
            physical_support='Qualitative compressibility effect only for extreme scenarios; conditional closure-based mechanics for controls',
            unresolved_assumptions='liquid EOS; shocks; shape stability; gas uniformity; fixed kappa and vapor; no mass transfer; finite molecular volume; continuum',
            negative_far_field_pressure=float(r['min_external_pressure_Pa'])<0,
            uncertainty='Experimental attainability and exact physical extrema are not established',
            next_stage_use=use,conditions=reason,
            category_note='Operational project admission categories, not universal physical cutoffs; no later stage run'))
    write_csv(RESULTS/'admission.csv',admitted)
    km=[r for r in rows if r['model']=='KM' and r['completed']=='True']
    rp=[r for r in rows if r['model']=='RP' and r['completed']=='True']
    extreme=max(km,key=lambda r:float(r['max_abs_U_m_s']))
    extreme_rp=max(rp,key=lambda r:float(r['max_abs_U_m_s']))
    incomplete=[r for r in rows if r['completed']=='False']
    c=read_csv(RESULTS/'c_limit.csv')+read_csv(RESULTS/'c_limit_extended.csv')
    all_errors={}
    for model in ('RP','KM'):
        for metric in ('Rmin_m','Rmax_m','max_abs_U_m_s','collapse_time_s','collapse_duration_s'):
            all_errors[model+'_'+metric]=max(float(r['relative_error_'+metric]) for r in cv
                                            if r['model']==model and float(r['rtol'])<=1e-10)
    # Locate the steepest finite-grid KM interval without declaring criticality.
    kr=sorted((r for r in km if float(r['Pa_bar'])>=1),key=lambda r:float(r['Pa_bar']))
    slopes=[]
    for a,b in zip(kr[:-1],kr[1:]):
        d=float(b['Pa_bar'])-float(a['Pa_bar'])
        slopes.append(dict(left_bar=float(a['Pa_bar']),right_bar=float(b['Pa_bar']),
            log_speed_slope_per_bar=float(np.log(float(b['max_abs_U_m_s'])/float(a['max_abs_U_m_s']))/d),
            log_radius_slope_per_bar=float(np.log(float(b['Rmin_m'])/float(a['Rmin_m']))/d)))
    write_csv(RESULTS/'scan_slopes.csv',slopes)
    steepest=max(slopes,key=lambda r:r['log_speed_slope_per_bar'])
    save_json(RESULTS/'audit_summary.json',dict(amplitudes=len(km),runs=len(rows),
        KM_extreme=extreme,RP_extreme_complete=extreme_rp,
        incomplete_runs=[{k:r[k] for k in ('model','Pa_bar','last_time_s','last_R_m','last_U_m_s','solver_message')} for r in incomplete],
        convergence_max_relative_errors=all_errors,steepest_KM_mesh_interval=steepest,
        final_c_limit=c[-1],physical_conclusion='Strongest KM result is a formal converged scenario, not an established physically attainable extremum'))
    tests=RESULTS/'test_results.txt'
    test_text=tests.read_text().strip() if tests.exists() else 'Tests are run after generation; see test_results.txt in the final snapshot.'
    env=json.loads((RESULTS/'clean_environment.json').read_text())
    text=f'''# Numerical and physical validation

Reference: accepted RP commit d59bc0d0119b514152f763e2537a0477c7c63663. The model file SHA256 is checked by an automated test. Only the two authorized MR1 documentation files change. The old MR1 checksum list remains the historical snapshot.

## Coverage and source checks

All five named primary papers and Part B of Plasma-plan.tex were read. Equation convention and full-source locations are documented in source_audit.md. The printed minus prefactor in Sadighi-Bonabi is not combined with the adopted plus-prefactor KM. Source hashes, solver/library versions and settings are saved in run_metadata.json.

{len(km)} amplitudes, {len(rows)} total RP/KM runs, {len(incomplete)} incomplete runs. All {len(km)} KM runs complete. No Mach filter, pressure clipping, velocity clipping or denominator regularization is used. All incomplete records retain their interval, state and reason; their full-interval comparison differences are undefined and their plotted lines have gaps.

## Automated and independent checks

Mandatory tests cover static gas balance, forcing phase, dimensional base-unit rescaling, the viscous pressure derivative along an independent quadratic path, the implicit residual, nonzero KM initial acceleration, algebraic RP limit, unchanged reference hash, positive radius/domain guards, event zero crossings, table/trajectory consistency, saved-parameter reproduction, censoring and cycle correspondence. The accepted MR1 tests also pass independently (14 tests).

Independent clean environment: Python {env['python']}; package versions and differing numerical settings are in clean_environment.json. All six original amplitudes are repeated for both equations. Separate four-configuration convergence uses DOP853/Radau at 1.0, 1.10, 1.15 and 1.30 bar. Convergence bounds are empirical observed differences, not rigorous error certificates. Other scan points are labeled production-only rather than independently converged.

Largest production/tight/Radau relative differences over these four anchors:

| Model | Rmin | Rmax | max abs(U) | collapse time | paired duration |
|---|---:|---:|---:|---:|---:|
'''
    for model in ('RP','KM'):
        vals=[all_errors[model+'_'+x] for x in ('Rmin_m','Rmax_m','max_abs_U_m_s','collapse_time_s','collapse_duration_s')]
        text+='| '+model+' | '+' | '.join(f'{x:.6g}' for x in vals)+' |\n'
    text+=f'''
The RP 1.10 bar production run has greater sensitivity than the other anchors; this is retained explicitly, with a numerical agreement criterion of 1e-4 relative across the checked metrics. KM uses 2e-5; its actual differences above are much smaller. Root-time tolerance is not physical timing accuracy. The strongest complete RP scan case receives an additional independent audit in strongest_RP_validation.csv.

## Mathematical incompressible limit

The exact explicit expression reduces to the imported RP right-hand side as c tends to infinity. The moderate five-cycle trajectory converges rapidly. The extreme RP core requires much larger c because tiny phase errors amplify later rebounds. At the last multiplier {float(c[-1]['c_factor']):.6g}, the event-aware maximum radius-trajectory error/R0 is {float(c[-1]['max_R_difference_over_R0']):.6g}, and speed-trajectory error/reference peak is {float(c[-1]['max_U_difference_over_reference_max']):.6g}. This is an observed numerical-resolution limitation of the sharp RP trajectory, not a physical sound speed. Event extrema approach RP much earlier. All factors, nonmonotonic intermediate errors and time shifts are retained; no convergence curve is silently trimmed.

## Extreme KM and physical status

Strongest saved KM: {float(extreme['Pa_bar']):g} bar, Rmin={float(extreme['Rmin_m'])*1e6:.9g} um, peak inward speed={float(extreme['Umax_inward_m_s']):.9g} m/s, Mach={float(extreme['max_wall_Mach']):.9g}. D/R remains positive; the minimum is {float(extreme['min_D_over_R']):.9g}. The pressure/bulk-modulus proxy is {float(extreme['max_liquid_pressure_change_over_bulk_modulus']):.9g}; this is not a computed density change.

Steepest measured KM log-speed interval is [{steepest['left_bar']:.3f}, {steepest['right_bar']:.3f}] bar on the selected grid. This does not establish a bifurcation or a unique threshold. Global extrema may change cycle; the cycle_comparison file is required for event delays.

Numerical convergence does not establish physical validity. Shape instability, gas-pressure uniformity, thermal/mass transfer, constant kappa, vapor pressure, excluded volume, continuum and EOS validity remain unresolved. Tensile liquid pressure does not prove impossibility; cavitation nucleation statistics are absent. The strongest KM extrema require a validated Gilmore/EOS comparison and gas-closure checks. The follow-up cases are 1.0, 1.15 and 1.3 bar; no unverified Gilmore model or CFD is included. No temperature or reaction is evaluated.

## Tests in final snapshot

```text
{test_text}
```

Every report page and all 13 PDF/PNG figure pairs are visually checked before publication. Publication and clean-snapshot verification are recorded separately; file availability and remote SHA are verified after pushing.
'''
    (RESULTS/'validation.md').write_text(text)


if __name__=='__main__': main()
