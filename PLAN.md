# PLAN — Numerical Renormalization of the Duffing–Van der Pol System

> Documento de trabajo. Objetivo final: un **paper en LaTeX (estilo APS/RevTeX, en inglés)** compilado a PDF,
> respaldado por un **repositorio GitHub ordenado** (`Hammonne/Duffing-Vanderpol`), con código reutilizable,
> notebooks limpios, figuras reproducibles y un **README en español**.
>
> Claude **no** hace `git commit` ni `git push`: sólo entrega los comandos para que Jimena los ejecute.
> Este archivo se va marcando con checklists `[ ]` → `[x]` a medida que se avanza.

---

## 0. Decisiones acordadas (congeladas)

- [x] **Ecuación combinada DVdP (canónica):**
  `dφ/dt = Π`, `dΠ/dt = -ω²φ - ε φ³ + α (1 - φ²) Π`
  con **ε (Duffing)** y **α (VdP) independientes**. Esto permite renormalizar cada efecto por separado
  (ε→0 recupera VdP puro; α→0 recupera Duffing puro) y luego combinarlos.
- [x] **Idioma:** paper en **inglés** (RevTeX/APS); README en **español**.
- [x] **Método numérico de la solución "full":** se hace **benchmark** (RK4 manual vs `solve_ivp`:
  RK45 / DOP853 / Radau) y se adopta el mejor para los tres sistemas, documentándolo.
- [x] **Renormalización UNIVERSAL (requisito clave):** un único método de renormalización debe resolver
  **Duffing solo, Van der Pol solo y DVdP combinado**, y aproximar bien la solución real **para cualquier
  parámetro, tanto suave (ε,α≈0.1) como fuerte (ε,α≈1 o mayores)**. No se asume régimen infinitesimal fijo.
- [x] **Casos de estudio:** SÍ se incluyen (Sección III del AFOSR) — benchmarking de la renormalización
  contra otros métodos y aplicación a control no lineal.
- [x] **Estructura del repo:** módulos `.py` reutilizables + notebooks limpios (uno por caso) + `paper/` + `figures/`.
- [x] **Erratum detectado:** el PDF escribe VdP como `ε(1+φ²)Π` (Ec. 17b); el estándar y los notebooks usan
  `ε(1-φ²)Π`. Se usa `(1-φ²)` y se deja nota en el paper.

### Marco teórico (de `project_UTEC_AFOSR.pdf`, Sección II)

- Oscilador armónico (orden ε⁰): `φ₀ = A cos(ϕ₀)`, `Π₀ = -ωA sin(ϕ₀)`, `ϕ₀ = ωt`.
- **Duffing (Caso 1):** la corrección de 1er orden `(φ¹,Π¹)` crece secularmente; la renormalización
  proyecta esa divergencia sobre un **corrimiento de frecuencia Δω** (generadores α_ω, β_ω).
  Para ω=1, ε: el shift de frecuencia `α_ω = 3/8` coincide con Bender–Orszag (multiple scales).
- **Van der Pol (Caso 2):** aquí NO basta cambiar frecuencia; se parametriza con **amplitud variable A(t)**
  (generador β_A). `dA_R/dt = ε β_A(A_R)`.
- **Combinado (nuevo, el aporte):** requiere **ambos** flujos simultáneos — Δω (de ε, Duffing) y ΔA
  (de α, VdP). La hipótesis de por qué "no funciona" actualmente es que los notebooks combinados usan
  un solo generador (sólo β_A o mezclan ε/α), por lo que no capturan los dos sectores a la vez.
- **Cost function / fitting:** los generadores se obtienen minimizando
  `I = Σ_A ∫ |φ¹_⊥,A|² dt`, i.e. proyectando `φ¹_A` sobre la base `∂φ₀_A/∂ϕ₀`, `∂φ₀_A/∂ω`, `∂φ₀_A/∂A`
  (forma cuadrática ⇒ mínimos cuadrados lineales).

### Referencias

- `[AFOSR]` project_UTEC_AFOSR.pdf — documento guía (Secciones II.A, II.B, II.C).
- `[GS21]` Gálvez Ghersi & Stein, Phys. Rev. E 104, 034219 (2021), arXiv:2106.08410 — algoritmo general ("paper de Thomas").
- `[Roy21]` Roy, Das, Banerjee, Int. J. Non-Linear Mech. 135 (2021) 103771 — referencia de método numérico (VdP–Mathieu–Duffing).
- Bender & Orszag (1999) — multiple scales, chequeo α_ω = 3/8.

---

## Inventario actual (notebooks del root)

| Archivo | Contenido | Destino |
|---|---|---|
| `Copy_of_Duffing_Renorm.ipynb` | Renorm. Duffing 1er orden + cuadraturas (trapz/Simpson/Gauss) + errores RMS + fases | Base del **notebook Duffing** limpio |
| `Duffing.ipynb` | RK4 manual + Duffing forzado/caótico + Lyapunov | Material auxiliar (apéndice / descartar forzamiento) |
| `Duffing_pert.ipynb` | Perturbación regular vs exacta (libro) + secular growth | Base de la sección **perturbativa** |
| `Duffing_Van_der_Pol___combinado.ipynb` | 4 variantes RG (VdP, Duffing, VdP, selector) — desordenado | Fuente para **combinado** |
| `Renato_DVdP.ipynb` | RG con `system_type` selector (duffing/vdp) + DVdP combinado | **Mejor base** del motor RG unificado |
| `Renormalizacion_VanDerPol_cambiado_para_duffing.ipynb` (x2, duplicado) | RG VdP adaptado a Duffing + extracción Ω(A) | Fuente VdP; **eliminar duplicado** |

---

## FASE 1 — Auditoría y consolidación del código (sin tocar física todavía) ✅

- [x] 1.1 Instalar Jupyter/nbconvert (`python -m jupyter`, `python -m nbconvert`).
- [x] 1.2 Ejecutar notebooks clave y registrar estado → `docs/audit.md`.
- [x] 1.3 Diagnosticado el bug del "merge" (Bugs A–D en `docs/audit.md`): un solo ε, generadores no combinados, perturbación mal planteada, paso RK4.
- [x] 1.4 Eliminado el duplicado exacto; originales respaldados en `legacy/`.

## FASE 2 — Esqueleto del repositorio ✅ (hecho: src/, notebooks/, paper/, figures/, docs/references/, tests/, legacy/, requirements.txt, .gitignore)

- [ ] 2.1 Crear estructura dentro de `Duffing-Vanderpol/`:
  ```
  src/dvdp/
    __init__.py
    systems.py        # define full systems: duffing, vdp, dvdp (ε, α independientes)
    integrators.py    # rk4 manual + wrappers solve_ivp (RK45/DOP853/Radau)
    background.py      # φ0, Π0 y derivadas ∂φ0/∂ϕ0, ∂φ0/∂ω, ∂φ0/∂A
    perturbation.py    # sistema de 1er orden (fuente = término de deformación)
    renormalization.py # fitting de generadores (β_ω, β_A), flujo de parámetros, φ_R
    quadrature.py      # trapz / Simpson / Gauss-Legendre (de Copy_of_Duffing_Renorm)
    metrics.py         # errores RMS, norma del error φ⊥
  notebooks/
    01_perturbation_vs_exact.ipynb
    02_duffing_renormalization.ipynb
    03_vanderpol_renormalization.ipynb
    04_combined_dvdp.ipynb
    05_numerical_method_benchmark.ipynb
  paper/   figures/   docs/   tests/
  ```
- [ ] 2.2 `requirements.txt` (numpy, scipy, matplotlib, jupyter, nbconvert) + `.gitignore` (`__pycache__`, `*.aux`, `*.log`, etc.).
- [ ] 2.3 Mover PDFs de referencia a `docs/references/`.

## FASE 3 — Núcleo de renormalización unificado (`src/dvdp`) ✅

- [x] 3.1 `systems.py`: `full_rhs(ε,α)`, `duffing_rhs`, `vdp_rhs`, `dvdp_rhs`.
- [x] 3.2 `background.py` + `perturbation.py`: fuentes separables `source_duffing=-φ0³`,
  `source_vdp=(1-φ0²)Π0`; superposición lineal. `quadrature.py`, `integrators.py` (RK4 corregido), `metrics.py`.
- [x] 3.3 `renormalization.py`: motor universal. Proyección robusta por **pendiente secular**
  (descomposición puntual en direcciones del background + ajuste lineal), generadores
  `β_A`, `β_θ` por deformación, combinación `β=ε·β^D+α·β^V`, flujo ODE no lineal, reconstrucción `φ_R,Π_R`.
- [x] 3.4 `tests/test_physics.py`: TODOS PASAN. Validado:
  - Duffing ω=1: `β_θ/A² = 0.375` (= **3/8**, Bender–Orszag) y `β_A≈0`.
  - VdP: `β_A = ½A(1-A²/4)` (ciclo límite A→2), `β_θ≈0`.
  - RG ≫ naive en Duffing/VdP/combinado, suave y fuerte (ε,α hasta 1).

**Resultados cuantitativos (RMS vs full, t∈[0,100]):**

| caso | ε | α | RG | naive |
|---|---|---|---|---|
| Duffing | 0.1 | – | 0.014 | 1.67 |
| Duffing | 1.0 | – | 1.07 | 15.5 |
| VdP | – | 0.1 | 0.030 | – |
| VdP | – | 1.0 | 2.09 | – |
| DVdP | 0.1 | 0.1 | 0.79 | 3.11 |
| DVdP | 1.0 | 1.0 | 1.91 | 22.3 |

## FASE 4 — Parte 1: reproducir el paper (Duffing + VdP independientes) ✅
- [x] 4.1 NB 01 + `fig_secular_growth.pdf`.
- [x] 4.2 NB 02 Duffing + `fig_duffing_compare.pdf` (RMS RG=0.014 ≪ naive=1.67).
- [x] 4.3 NB 03 VdP + `fig_vdp_compare.pdf` (ciclo límite A→2).
- [x] 4.4 Validado numéricamente: `β_θ/A²=3/8` (Duffing), `β_A=½A(1-A²/4)` (VdP) + `fig_generators.pdf`.

## FASE 5 — Benchmark del método numérico ✅
- [x] 5.1 `scripts/benchmark_integrators.py` → `docs/benchmark.md` + `fig_benchmark.pdf` (NB 05).
- [x] 5.2 Bug de paso del RK4 manual corregido en `integrators.rk4`.
- [x] 5.3 **DOP853** elegido (rtol=1e-10), fijado como default en `integrators.solve_full`, justificado en paper+README.

## FASE 6 — Ecuación combinada DVdP (el aporte nuevo) ✅
- [x] 6.1 NB 04 DVdP full (barrido ε, α).
- [x] 6.2 Renormalización combinada por superposición de generadores (`β=ε·β^D+α·β^V`); causa del fallo previo documentada en `docs/audit.md`.
- [x] 6.3 `fig_dvdp_soft.pdf`, `fig_dvdp_strong.pdf` (temporal + espacio de fases).
- [x] 6.4 Límites ε→0 ⇒ VdP, α→0 ⇒ Duffing verificados en tests.
- [x] 6.5 Rango de validez cuantificado en `fig_error_scaling.pdf`.

## FASE 7 — Paper LaTeX (RevTeX, inglés) ✅
- [x] 7.1 `paper/main.tex` (revtex4-2, aps/pre): I Intro · II Method (A/B/C) · III Combined DVdP · IV Benchmark · V Cases of study · VI Discussion.
- [x] 7.2 Ecuaciones transcritas y limpias.
- [x] 7.3 Figuras insertadas con captions APS.
- [x] 7.4 `paper/refs.bib` (12 refs AFOSR + Roy + SciPy).
- [x] 7.5 Compilado: `main.pdf` (5 págs), 0 errores, sin citas/refs sin resolver.

## FASE 8 — README y documentación ✅
- [x] 8.1 README en español (motivación, ecuaciones, estructura, reproducción, método, limitaciones).
- [x] 8.2 `docs/audit.md`, `docs/benchmark.md`, `docs/references/`.

## FASE 9 — Entrega
- [x] 9.1 Verificación: tests pasan, notebooks ejecutan, `main.pdf` compila.
- [ ] 9.2 Entregar a Jimena la secuencia de comandos git (ver final del chat).

---

## Preguntas abiertas / riesgos  (RESUELTAS)

- [x] Jupyter → **instalar** (autorizado). Se instala jupyter/nbconvert/nbclient/ipykernel.
- [x] Parámetros: el método debe funcionar para **cualquier** ε, α (suave y fuerte); figuras del combinado
  incluyen barrido suave (0.1) **y** fuerte (≥1). Riesgo técnico: el ansatz lineal infinitesimal del AFOSR
  se rompe para ε,α grandes ⇒ se implementan **generadores no lineales** β_A(A), β_ω(ω) resueltos como ODE
  del flujo de parámetros (no linealizados), que es lo que permite capturar el régimen fuerte.
- [x] El paper **SÍ incluye** Sección III (casos de estudio y benchmarking) y IV (discusión).

## Nota de diseño: por qué "universal" es factible

El método de `[GS21]` (Thomas) es general para sistemas con N grados de libertad: el núcleo no depende
de la forma del término de deformación, sólo de (a) la base de derivadas del background
`{∂φ0/∂ϕ0, ∂φ0/∂ω, ∂φ0/∂A}` y (b) el término fuente de 1er orden. Por eso un mismo `renormalization.py`
puede recibir `source = -ε φ0³` (Duffing), `source = α(1-φ0²)Π0` (VdP) o la suma (DVdP) y, activando los
generadores pertinentes (Δω y/o ΔA), produce la solución renormalizada en cada caso. Para régimen fuerte se
resuelve el flujo completo no lineal β(λ) en lugar de truncar a ∆t lineal.
