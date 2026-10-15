# Auditoría de notebooks (Fase 1)

Fecha: 2026-10-15. Entorno: Python 3.10, numpy/scipy/matplotlib OK, Jupyter/nbconvert instalados
(`python -m jupyter`, `python -m nbconvert`). Backend de ejecución: `MPLBACKEND=Agg`.

## Estado de ejecución

| Notebook | ¿Corre? | Notas |
|---|---|---|
| `Renato_DVdP.ipynb` | ✅ Sí | Motor RG más completo (A y θ como parámetros, dos generadores β_A, β_θ). |
| `Copy_of_Duffing_Renorm.ipynb` | ✅ Sí | Implementación más fiel al paper: cost-function con matriz normal M y cuadraturas (trapz/Simpson/Gauss). |
| `Renormalizacion_VanDerPol_cambiado_para_duffing.ipynb` | ⚠️ Corre pero con bug conceptual | Ver "Bug A". Duplicado exacto: `...(1).ipynb` → eliminar. |
| `Duffing_pert.ipynb` | ✅ (perturbación regular) | Material de secular growth; contiene bug de paso RK4 (Bug C). |
| `Duffing.ipynb` | ✅ | Duffing forzado/caótico + Lyapunov; fuera del alcance del core (apéndice). |
| `Duffing_Van_der_Pol___combinado.ipynb` | Parcial | 4 variantes desordenadas, celdas vacías; superado por `Renato_DVdP`. |

## Motores de renormalización encontrados (2 estilos)

1. **Estilo "cost-function / matriz normal"** (`Copy_of_Duffing_Renorm`): proyecta la perturbación
   `(φ¹,Π¹)` sobre la base de derivadas del background `{∂φ0/∂ϕ0, ∂φ0/∂ω}` resolviendo el sistema
   normal `M Δλ = b`, con `M,b` calculados por cuadratura (trapz/Simpson/Gauss). **Es la implementación
   fiel de la Ec. (16) del AFOSR** (minimización de `I = Σ_A ∫|φ¹_⊥|²dt`). → base para `renormalization.py`.

2. **Estilo "ajuste de pendiente secular"** (`Renato_DVdP`, `Renorm_VanDerPol`): ajusta β como la pendiente
   del crecimiento secular, proyectando `φ¹` sobre `{t·cosωt, -A·t·sinωt}` por mínimos cuadrados.
   Más simple; captura bien la tendencia secular. `Renato_DVdP` generaliza con A y θ.

Ambos son equivalentes al orden dominante; adoptaremos (1) como núcleo (más riguroso y da Δω además de ΔA),
validándolo contra (2).

## CAUSA RAÍZ del "no funciona" al combinar Duffing + VdP

### Bug A — un solo parámetro de deformación en vez de ε y α independientes
En `Renato_DVdP` (celda DVdP) y `Renorm_VanDerPol`, el término fuente combinado se escribe
`source = -p0³ + (1-p0²)·pi0` y luego el flujo usa **un único** `epsilon` (`dA=ε·β_A`, `dθ=ε·β_θ`).
No hay forma de encender Duffing y VdP con intensidades distintas. El usuario quiere ε (Duffing) y
α (VdP) **independientes**. Esto es el bloqueador principal.

### Bug B — no se combinan los dos sectores de generadores
- Duffing renormaliza **frecuencia/fase** (Δω, β_θ).
- VdP renormaliza **amplitud** (ΔA, β_A).
- El combinado necesita **ambos a la vez**. `Renorm_VanDerPol` sólo extrae β_A (ignora la fase) y lo
  aplica al sistema DVdP → captura la amplitud pero no el corrimiento de frecuencia del término cúbico.
  `Renato_DVdP` sí usa ambos, pero atados a un solo ε (Bug A).

### Bug C — sistema perturbativo mal planteado en `Renorm_VanDerPol`
Su `perturbative_system` reintroduce la **no linealidad completa**
(`dPi1 = -ω²φ1 - εφ1³ + α(1-φ1²)Pi1`) en vez de la ecuación **lineal de primer orden** con fuente
construida desde el background `p0,pi0`. El primer orden correcto (Ec. 6b / 19b del AFOSR) es
`dΠ1/dt = -ω²φ1 + source(p0,pi0)`. Con IC cero, su versión da φ1≈0 (fuente efectiva nula) ⇒ resultados
"planos"/no favorables. `Renato_DVdP` y `Copy_of_Duffing_Renorm` sí lo plantean bien.

### Bug D — paso temporal mal en RK4 manual (`Duffing.ipynb`, `Duffing_pert.ipynb`)
`h=(tf-t0)/len(t_eval)` con `t_span=(0,300)` pero `t_eval=linspace(0,3000,...)`: desalineación t_span/t_eval
y paso inconsistente. Descartaremos el RK4 manual salvo para el benchmark (Fase 5), donde se corrige.

## SOLUCIÓN propuesta (método universal) — para Fase 3

Como la ecuación de primer orden es **lineal**, vale superposición:
`φ¹[ε·Duffing + α·VdP] = ε·φ¹_D + α·φ¹_V`.
Por tanto se calculan **generadores separados** proyectando cada fuente pura:
- `β^D(λ)` desde `source_D = -(φ0)³`  (contribuye principalmente a β_ω)
- `β^V(λ)` desde `source_V = (1-φ0²)·Π0`  (contribuye principalmente a β_A)
y el flujo combinado es
`dλ/dt = β^(0)(λ) + ε·β^D(λ) + α·β^V(λ)`, con `λ=(A, θ/ω)`.

Esto da un **único método universal**: ε=0 ⇒ VdP puro; α=0 ⇒ Duffing puro; ambos ⇒ DVdP.
Para **régimen fuerte** (ε,α grandes) se resuelve el **flujo no lineal completo** `dλ/dt=…β(λ)` (β como
función interpolada de λ sobre una grilla), en lugar de truncar a `Δt` lineal — esto captura la parte
no infinitesimal, que es la limitación que el propio AFOSR señala.

## Acciones de limpieza

- [x] Confirmado duplicado exacto `Renormalizacion_VanDerPol_cambiado_para_duffing (1).ipynb`.
- [ ] Eliminar duplicado.
- [ ] Mover PDFs a `docs/references/`.
