# Convergencia electrónica separada de los otros errores

**Las tres mallas electrónicas pasan el control de convergencia al separar
el regulador y la interpolación fonónica fija.** No basta con observar que el
error combinado está bajo tolerancia: se completaron las referencias BCS con
η=10⁻⁸ para 1260 y 2520 estados y se compararon los operadores ya guardados.
El único cálculo nuevo tardó **4,26 s en Geminga**, con un hilo.

| Estados | Error combinado contra la referencia ideal | Error electrónico separado |
|---:|---:|---:|
| 630 | 2,99371×10⁻⁴ | 2,91333×10⁻⁴ |
| 1260 | 2,82810×10⁻⁴ | 7,28419×10⁻⁵ |
| 2520 | 2,91511×10⁻⁴ | 2,18741×10⁻⁵ |

La segunda columna conserva el resultado histórico: **no es una serie
monótona y no se declara como tal**. La tercera mantiene las ocupaciones
electrónicas nativas y sus actividades de Pauli; sólo evalúa Bose directamente
en la energía de cada transición y emplea la referencia BCS con la misma η.
Es la comparación en que la cuadratura electrónica se refina por separado,
como pide el contrato.

Las **18 series** —tres campos, tres perfiles y dos canales— disminuyen en
ambos refinamientos. Ningún caso requiere una excepción por redondeo ni una
tolerancia nueva. El mayor presupuesto numérico de las referencias es
2,53524×10⁻⁷, muy por debajo de su límite original 2×10⁻⁴.

La separación se hizo sobre los vectores completos del RHS. La diferencia
total es exactamente la suma de las diferencias por interpolación Bose,
reconstrucción Pauli, integración electrónica y regulador BCS. El residuo
normalizado de esta identidad fue como máximo 4,05×10⁻²¹. No se sumaron normas
para fingir que los errores se cancelan ni se corrigió el operador.

Esta evidencia permite cerrar el control **electrónico separado de la matriz
histórica Ωmín=0,01**. La comparación estática seleccionada Ωmín=0,005,
630/1025 ya pasa en `selected_continuum.json`, pero estos archivos no crean
tres nuevas mallas a ese corte ni establecen un dictamen temporal o global.

Archivos nuevos:

- `electronic_three_grid_decomposition.json`: todas las series, componentes,
  presupuestos de referencia y comprobaciones de disminución.
- `same_eta_bcs_reference_2.json` y `same_eta_bcs_reference_4.json`: las dos
  referencias faltantes, con integrales independientes de orden 4/32 y 8/64.
- `sandbox/stage2_cells/resume_20260921/audit_electronic_three_grids.py`:
  constructor reproducible de esta evidencia, sin crear redes ni dinámica.

Comando reproducible desde `/home/jdiaz/pysnspd`:

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  sandbox/stage2_cells/resume_20260921/audit_electronic_three_grids.py
```
