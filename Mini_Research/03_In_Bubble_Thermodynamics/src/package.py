"""Assemble uncertainty notes, validation summary and self-contained file manifests."""
from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import pandas as pd
from .bubble_thermo import ROOT,MR2,MECHANICAL_COMMIT
from .run_experiments import json_write

OUT=ROOT/'results'


def uncertainty():
    a=pd.read_csv(OUT/'adiabatic_summary.csv');b=pd.read_csv(OUT/'thermal_summary.csv')
    c=pd.read_csv(OUT/'solver_convergence.csv');q=pd.read_csv(OUT/'input_grid_convergence.csv')
    upstream=pd.read_csv(MR2/'results/convergence.csv');rows=[]
    sensitivity=pd.read_csv(OUT/'adiabatic_sensitivity.csv')
    a['analysis']='reference_mechanics'
    sensitivity['analysis']='frozen_gamma_T0_sensitivity'
    for r in pd.concat([a,sensitivity],ignore_index=True).to_dict('records'):
        z=upstream[(upstream.model==r['mechanics'])&(upstream.Pa_bar==r['Pa_bar'])&
                   (upstream.method=='DOP853')&(upstream.rtol==1e-10)]
        error=(1-float(z.iloc[0].relative_error_Rmin_m))**(-3*(r['gamma']-1))-1 if len(z) else np.nan
        rows.append(dict(mechanics=r['mechanics'],thermodynamics=r['thermodynamics'],Pa_bar=r['Pa_bar'],
            gas=r['gas'],gamma=r['gamma'],T0_K=r['T0_K'],analysis=r['analysis'],
            eta=0,observed_Tmax_K=r['observed_Tmax_K'],t_Tmax_s=r['t_Tmax_s'],
            completed=r['completed'],applicability=r['applicability'],
            numerical_peak_difference_relative=error,
            interpolation_peak_difference_relative=r['refined_vs_interpolant_relative_error'],
            numerical_note='Propagated upstream radius difference; not a rigorous bound' if len(z) else
                'No amplitude-specific tighter reference; inherited production-only status where so classified',
            physical_note=('Frozen mechanics and water; effective gamma is not a chemical identity; '
                'changed pressure is inconsistent with mechanical closure unless gamma=1.4' if
                r['analysis']=='frozen_gamma_T0_sensitivity' else
                'Unquantified closure uncertainty; no experimental temperature validation')))
    for r in b.to_dict('records'):
        df=c if r['thermodynamics']=='B_fully_coupled' else q
        z=df[(df.Pa_bar==r['Pa_bar'])&(df.gas==r['gas'])&(df.eta==r['eta'])]
        if r['thermodynamics']=='B_fully_coupled':z=z[z.configuration=='production']
        else:z=z[(z.grid=='full_export_production')&(z.thermodynamics==r['thermodynamics'])]
        grid=q[(q.Pa_bar==r['Pa_bar'])&(q.gas==r['gas'])&(q.eta==r['eta'])&
               (q.grid=='export_stride2')&(q.thermodynamics==r['thermodynamics'])]
        rows.append(dict(mechanics='KM',thermodynamics=r['thermodynamics'],Pa_bar=r['Pa_bar'],gas=r['gas'],
            eta=r['eta'],gamma=r['gamma'],T0_K=r['T0_K'],analysis='thermal_model',
            observed_Tmax_K=r['observed_Tmax_K'],t_Tmax_s=r['t_Tmax_s'],completed=r['completed'],
            applicability=r['applicability'],numerical_peak_difference_relative=float(z.iloc[0].relative_Tmax_error),
            interpolation_peak_difference_relative=float(grid.iloc[0].relative_Tmax_error) if len(grid) else np.nan,
            numerical_note='Production versus tighter tolerance and half maximum step; grid sensitivity separate',
            physical_note='Constant ambient conductivity and heat capacity; fixed vapor; uniform ideal gas; no validation'))
    pd.DataFrame(rows).to_csv(OUT/'numerical_uncertainty.csv',index=False)


def summary():
    a=pd.read_csv(OUT/'adiabatic_summary.csv');b=pd.read_csv(OUT/'thermal_summary.csv')
    c=pd.read_csv(OUT/'solver_convergence.csv');q=pd.read_csv(OUT/'input_grid_convergence.csv')
    r=pd.read_csv(OUT/'adiabatic_recovery.csv')
    v=dict(accepted_mechanical_commit=MECHANICAL_COMMIT,adiabatic_cases=len(a),thermal_cases=len(b),
        full_thermal_cases=int(b.completed.sum()),censored_mechanical_cases=int((~a.completed).sum()),
        max_adiabatic_refined_agreement=float(a.refined_vs_interpolant_relative_error.max()),
        max_coupled_production_peak_difference=float(c[c.configuration=='production'].relative_Tmax_error.max()),
        max_independent_Radau_peak_difference=float(c[c.configuration=='independent_Radau'].relative_Tmax_error.max()),
        all_current_convergence_runs_completed=bool(c.completed.all()),
        max_one_way_tolerance_peak_difference=float(q[(q.grid=='full_export_production')&
            (q.thermodynamics=='B_one_way_post_processing')].relative_Tmax_error.max()),
        max_fine_input_peak_difference=float(q[(q.grid=='export_stride2')&
            (q.thermodynamics=='B_one_way_post_processing')].relative_Tmax_error.max()),
        max_scaled_production_gas_energy_residual=float(b.energy_residual_scaled.max()),
        max_independent_quadrature_energy_residual=float(max(c.gauss_energy_scaled_residual.max(),q.gauss_energy_scaled_residual.max())),
        max_coupled_adiabatic_peak_difference_from_accepted=float(r.relative_Tmax_vs_saved_mechanics.max()),
        source_zero_for_all_runs=bool((b.source_energy_J==0).all()),
        physical_validation=False,acceptance='awaiting independent review')
    json_write(OUT/'validation_summary.json',v)
    text='''# Проверки Mini-Research 3

Исследование подготовлено для независимой проверки; окончательное принятие не заявляется.
Результаты тестов находятся в `tests.txt`, проверки чистого окружения в
`clean_environment_verification.json`, численные показатели в `validation_summary.json`.

'''
    text+='| Проверка | Результат |\n|---|---|\n'
    for key,value in v.items():text+=f'| {key} | {value} |\n'
    text+='''
Численные разности показывают устойчивость выбранных расчётов и не являются
строгими границами погрешности. Для one-way отдельно учитываются погрешности
энергетического интегрирования, новой интерполяции и унаследованной механики.
Сохранённые механические экстремумы добавлены к каждой прореженной сетке;
проверка длительности и отложенных узлов оценивает поведение между ними.

Все предшествующие файлы проверены по SHA256 и принятому коммиту. Неполные RP
результаты не превращены в глобальные за пять циклов. Завершённые экстремальные
RP результаты оставлены диагностическими. Источники, конвенция KM и физические
ограничения описаны в `source_audit.md` и научном отчёте.
'''
    (OUT/'validation.md').write_text(text)


def manifests():
    ignored={'__pycache__','.venv','.pytest_cache'}
    aux={'.aux','.log','.out','.fls','.fdb_latexmk','.gz'}
    files=[]
    for p in ROOT.rglob('*'):
        rel=p.relative_to(ROOT)
        if not p.is_file() or any(x in ignored for x in rel.parts):continue
        if p.name in {'FILE_LIST.txt','SHA256SUMS'}:continue
        if rel.parts[0]=='report' and p.suffix in aux:continue
        if p.suffix=='.log':continue
        files.append(p)
    files=sorted(files)
    names=[p.relative_to(ROOT).as_posix() for p in files]
    (ROOT/'FILE_LIST.txt').write_text('\n'.join(sorted(names+['FILE_LIST.txt','SHA256SUMS']))+'\n')
    files.append(ROOT/'FILE_LIST.txt')
    (ROOT/'SHA256SUMS').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(ROOT).as_posix()+'\n' for p in sorted(files)))


def main():
    uncertainty();summary();manifests()


if __name__=='__main__':main()
