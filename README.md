# Duffing–Van der Pol: Renormalización numérica universal

Un método de **renormalización numérica (DRG)** que aproxima la solución del
oscilador armónico deformado por un término de **Duffing** (rigidez cúbica) y/o
un término de **Van der Pol** (amortiguamiento no lineal), usando únicamente la
respuesta perturbativa de primer orden.

La ecuación maestra (con dos parámetros de deformación **independientes**
ε = Duffing y α = Van der Pol) es

```
dφ/dt = Π
dΠ/dt = -ω²·φ - ε·φ³ + α·(1 - φ²)·Π
```

- `ε = 0` → Van der Pol puro
- `α = 0` → Duffing puro
- ambos ≠ 0 → sistema combinado **DVdP**

Un **único motor** (`src/dvdp`) resuelve los tres casos: ese es el aporte
central del trabajo.

## Idea del método

1. El término no lineal se trata como una deformación de orden ε/α. La
   corrección de primer orden crece **secularmente** (diverge en el tiempo).
2. La renormalización reabsorbe ese crecimiento promoviendo la **amplitud `A`**
   y la **fase `θ`** del fondo `φ⁰ = A·cos(ωt+θ)` a funciones lentas del tiempo.
3. Los **generadores** `β_A(A)` y `β_θ(A)` se obtienen proyectando el
   crecimiento secular sobre las direcciones tangentes del fondo (extracción
   robusta por **pendiente secular**).
4. Como la ecuación de primer orden es **lineal**, los generadores se
   **superponen**: `β = ε·β^Duffing + α·β^VdP`. Por eso el método es universal.
5. Se integra el flujo `dA/dt = β_A`, `dθ/dt = β_θ` (ODE no lineal completa, lo
   que permite acoplamiento fuerte) y se reconstruye `φ_R = A(t)·cos(ωt+θ(t))`.

### Resultados físicos validados
- **Duffing** (ω=1): corrimiento de frecuencia `β_θ = (3/8)·A²` — coincide con
  el resultado de multiple-scales de Bender–Orszag. La amplitud se conserva.
- **Van der Pol**: flujo de amplitud `β_A = ½·A·(1 - A²/4)` — ciclo límite en
  `A = 2`. Sin corrimiento de frecuencia.
- **Combinado**: ambos efectos a la vez; `ε=0` y `α=0` reducen a los casos puros.

## Estructura del repositorio

```
src/dvdp/              # motor (fuente única de la física)
  systems.py           # RHS full: Duffing / VdP / DVdP
  background.py        # fondo armónico y sus derivadas
  perturbation.py      # ecuación lineal de 1er orden + fuentes separables
  renormalization.py   # motor universal: generadores + flujo + reconstrucción
  integrators.py       # RK4 (corregido) y wrappers de solve_ivp; solución full/naive
  quadrature.py        # trapz / Simpson / Gauss-Legendre
  metrics.py           # errores RMS
notebooks/             # un notebook limpio por caso (01–05)
scripts/               # make_figures.py, benchmark_integrators.py, build_notebooks.py
paper/                 # artículo en LaTeX (RevTeX/APS) → main.pdf
figures/               # figuras generadas (PDF)
docs/                  # audit.md, benchmark.md, references/ (PDFs)
tests/                 # test_physics.py (regresiones físicas)
legacy/                # notebooks originales (respaldo)
```

## Cómo reproducir

```bash
pip install -r requirements.txt

# 1. Tests de física (deben pasar todos)
python tests/test_physics.py

# 2. Benchmark de integradores → docs/benchmark.md + figures/fig_benchmark.pdf
python scripts/benchmark_integrators.py

# 3. Todas las figuras del paper → figures/*.pdf
python scripts/make_figures.py

# 4. (Re)generar y ejecutar los notebooks
python scripts/build_notebooks.py
jupyter nbconvert --to notebook --execute --inplace notebooks/*.ipynb

# 5. Compilar el paper
cd paper
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

## Método numérico

El benchmark (ver `docs/benchmark.md`) compara RK4 de paso fijo contra RK45,
DOP853 y Radau de SciPy sobre los tres sistemas. Se adopta **DOP853**
(rtol=1e-10) como integrador de referencia por su mejor precisión por
evaluación en los casos no rígidos; Radau queda como respaldo si aparece rigidez
para α muy grande.

## Limitaciones

Para acoplamiento **muy fuerte** (ε, α ≳ 1) la aproximación captura la
envolvente de amplitud y el tamaño del ciclo límite, pero aparece una deriva de
fase a tiempos largos: el ansatz de un solo armónico no representa la forma
anarmónica del ciclo límite real. Es la limitación esperada de un esquema a
primer orden.

## Referencias

- `docs/references/project_UTEC_AFOSR.pdf` — documento guía (Canahuire, Gálvez
  Ghersi, Sánchez), Sección II.
- Gálvez Ghersi & Stein, Phys. Rev. E **104**, 034219 (2021).
- Roy, Das, Banerjee, Int. J. Non-Linear Mech. **135**, 103771 (2021).
- Bender & Orszag, *Advanced Mathematical Methods* (1999).
