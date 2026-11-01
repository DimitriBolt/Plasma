# MR8 independent scientific and editorial review

Marker: **MR8-REVIEW-083fc4a4fc635d0191b4c78b6b96a7ecf1e6b4f2**

Scientific publication SHA: `083fc4a4fc635d0191b4c78b6b96a7ecf1e6b4f2`.
Remote `main` was independently read back at this SHA after publication.
Both remotely downloaded PDFs match the locally verified SHA-256 fingerprints.
Subsequent receipt commits change only supporting publication/delivery records.

- [Scientific commit](https://github.com/DimitriBolt/Plasma/commit/083fc4a4fc635d0191b4c78b6b96a7ecf1e6b4f2)
- [Mini-Research 8 Technical Report, 9 pages](https://github.com/DimitriBolt/Plasma/blob/083fc4a4fc635d0191b4c78b6b96a7ecf1e6b4f2/Mini_Research/08_Physical_Feasibility/report/08_Physical_Feasibility.pdf)
- [Consolidated Final Scientific Report, 24 pages](https://github.com/DimitriBolt/Plasma/blob/083fc4a4fc635d0191b4c78b6b96a7ecf1e6b4f2/Final_Report/Consolidated_Final_Scientific_Report.pdf)
- [Source, data and reproduction instructions](https://github.com/DimitriBolt/Plasma/tree/083fc4a4fc635d0191b4c78b6b96a7ecf1e6b4f2/Mini_Research/08_Physical_Feasibility)
- [Validation journal](https://github.com/DimitriBolt/Plasma/blob/083fc4a4fc635d0191b4c78b6b96a7ecf1e6b4f2/Mini_Research/08_Physical_Feasibility/results/validation.md)
- [Standalone-reader assessment](https://github.com/DimitriBolt/Plasma/blob/083fc4a4fc635d0191b4c78b6b96a7ecf1e6b4f2/Mini_Research/08_Physical_Feasibility/results/standalone_reader_test.md)
- [Primary-source audit](https://github.com/DimitriBolt/Plasma/blob/083fc4a4fc635d0191b4c78b6b96a7ecf1e6b4f2/Mini_Research/08_Physical_Feasibility/results/source_audit.md)

## New calculations

1. Eighteen bounded one-cycle mechanical integrations: KM and Gilmore–Tait at
   1.00, 1.15 and 1.30 bar, each with production/strict DOP853 and independent
   Radau. Common MR2 polytropic gas, reference state and forcing isolate the
   liquid-model change. Tait n=7.15; B=rho0*c0²/n−p0=306.524 MPa matches the
   existing ambient acoustic state. Full viscous acceleration derivative and
   enthalpy derivative are included. No temperature is transplanted.
2. Read-only gas/liquid/Au/summed thermal identities from selected public MR7
   references, plus independent nonlinear caloric delta-U from saved paired
   delta-T and its original work/heat ledgers.
3. Ideal instantaneous all-absorption gas heating at fixed volume and fixed
   composition, inverting the existing NASA caloric integral only in its domain;
   separately, transport times and the slip/Reynolds requirements of prescribed
   constant-gap retention. No new chemistry, ionization or fusion state is invented.

At 1.30 bar, one-cycle Rmin is 0.550809 µm (KM) versus 0.513361 µm (Gilmore),
a 6.80% decrease. Peak inward speeds are 1072.98 versus 1031.29 m/s. Gilmore's
maximum local Mach is 0.504 and Tait density ratio 1.436. These remain formal
mechanical sensitivity results, not validated physical extrema or new gas temperatures.

For one nominal particle, Eabs=34.3373 fJ. The ideal caloric allowances at the
principal wet baseline peaks are 0.03347 and 0.02499 K; they are explicitly not
dynamic/universal bounds. The unchanged MR7 integrated global increments are
1.408816 and 2.184879 µK. The 1.15 bar pair completes numerically but remains
physically formal; its late event-centered sensitivity is unresolved (89.6%
mesh change). Wet 1.30 bar remains coefficient-censored, with no global increment.
Dry formal MR7 peaks 6263.58 and 22808.7 K and much larger invalid RP diagnostics
remain visible with distinct model and interval status.

## Mathematical, numerical and artifact checks

- Twelve tests pass in the working package and in an exact clean Git archive of
  the published scientific commit. This includes all 1491 preexisting tracked
  files being byte-identical; no MR1–MR7 result, source, root plan or historical
  artifact is changed.
- New equilibrium, Tait enthalpy integral/derivative, full implicit viscous
  equation, RP limit, gas domain and signed energy tests pass.
- Maximum mechanical control differences: 4.27e-10 relative Rmin, 3.70e-9
  relative peak speed, 1.47e-16 s selected-event difference. These are numerical
  discrepancies, not physical accuracy or universal error bounds.
- Summed thermal residual at most 1.35e-6 of the thermal scale; independent
  paired gas residual below 4.61e-8 of absorption. Net gas heat includes feedback
  and is not a tagged photon-delivery fraction.
- The initial clean staged-tree archive passes 12 tests, 80 manifest entries,
  fresh strong-Gilmore replay and both LaTeX builds; 15 figure/table files
  regenerate byte-identically, and both rebuilt PDFs have matching extracted
  text. The exact published scientific commit additionally passes 12 tests and
  all 82 manifest entries (64 MR8 + 18 Final); the extra files are verification
  records. The fresh cycle is not represented as a replay of all MR7 sweeps.
- Final Report has 20 pages from introduction through conclusions plus title,
  contents, references and the seven-report guide. All 24+9 PDF pages inspected;
  fonts embedded; no overfull boxes or undefined references. All seven exact
  supporting GitHub PDF links return 200 and readable PDF headers.
- Full named Keller–Miksis, Sadighi-Bonabi/Gheshlaghi, Flannigan/Suslick, NASA
  and photothermal sources checked. Exact Toegel article read from the preserved
  cache when fresh download returned 403. Gilmore implementation is grounded in
  the complete Denner primary manuscript, not an inaccessible original reprint.

The Final Report explicitly states three goals, all essential equations and
reference parameters, baseline versus laser results, formal extrema, physical
limits, thermonuclear prerequisites and specific missing tests. Footnotes give
full MR titles, verified section names and clickable exact PDF filenames.
Hahn is the intended reader; Gabitov remains the academic supervisor.

## Main scientific conclusion and unresolved premises

The nominal single-particle path changes the computed peak very little, but no
universal efficiency or physically reliable upper temperature has been derived.
The program does not substantiate thermonuclear feasibility, including rare
extreme configurations; unresolved configurations are not declared impossible.
Uniform gas temperature is not Ti, Te or a measured hot-region maximum.
No isotope-resolved fuel density, hot volume, ion-state history or reacting
lifetime exists, so DD yield is left undetermined rather than assigned arbitrary
inputs or zero. Water, acid spectroscopy and deuterated-acetone literature are
separated explicitly.

Remaining physical issues: tensile stability and competing nuclei; compressible
liquid EOS, shocks/nonsphericity; reacting/nonideal and nonequilibrium gas;
vapor/phase enthalpy; finite-Au moving local field and temperature-dependent
optics; realizable particle retention; and fuel/ion density–volume–time data.

Reproducibility is complete for the bounded new work using the whole tracked
repository and pinned Python environment. The public compact MR7 references
suffice. Its approximately 1 GB full local archive and unpublished history were
not uploaded or made ancestors.

Please independently assess the new equations, numerical evidence, energy
interpretation, status classification and Final Report's standalone readability
for Prof. Dr. Hahn. MR8 is not declared finally accepted by Codex. Scientific
and editorial acceptance remains with Mini-Research 8 (Web); official closure
of the eight-stage program remains with RoadMap. No MR9 is initiated.
