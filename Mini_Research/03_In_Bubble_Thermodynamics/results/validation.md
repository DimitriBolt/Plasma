# Проверки Mini-Research 3

Исследование подготовлено для независимой проверки; окончательное принятие не заявляется.
Результаты тестов находятся в `tests.txt`, проверки чистого окружения в
`clean_environment_verification.json`, численные показатели в `validation_summary.json`.

| Проверка | Результат |
|---|---|
| accepted_mechanical_commit | f1f3e812d77fd7fa058ce111ff484b1c8658f937 |
| adiabatic_cases | 112 |
| thermal_cases | 30 |
| full_thermal_cases | 30 |
| censored_mechanical_cases | 7 |
| max_adiabatic_refined_agreement | 6.4213079298269796e-12 |
| max_coupled_production_peak_difference | 6.714637734717144e-10 |
| max_independent_Radau_peak_difference | 1.9514290183764158e-10 |
| all_current_convergence_runs_completed | True |
| max_one_way_tolerance_peak_difference | 1.436052832204382e-09 |
| max_fine_input_peak_difference | 8.543969953755949e-08 |
| max_scaled_production_gas_energy_residual | 6.541570947982134e-06 |
| max_independent_quadrature_energy_residual | 3.703448742920665e-07 |
| max_coupled_adiabatic_peak_difference_from_accepted | 1.1241718667065469e-10 |
| source_zero_for_all_runs | True |
| physical_validation | False |
| acceptance | awaiting independent review |

Численные разности показывают устойчивость выбранных расчётов и не являются
строгими границами погрешности. Для one-way отдельно учитываются погрешности
энергетического интегрирования, новой интерполяции и унаследованной механики.
Сохранённые механические экстремумы добавлены к каждой прореженной сетке;
проверка длительности и отложенных узлов оценивает поведение между ними.

Все предшествующие файлы проверены по SHA256 и принятому коммиту. Неполные RP
результаты не превращены в глобальные за пять циклов. Завершённые экстремальные
RP результаты оставлены диагностическими. Источники, конвенция KM и физические
ограничения описаны в `source_audit.md` и научном отчёте.
