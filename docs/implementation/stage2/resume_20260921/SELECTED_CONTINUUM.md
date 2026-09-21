# Comparación continua de la configuración seleccionada

**El control estático pasa con el corte real Ωmín=0,005.** Se compararon los
resultados ya guardados para 630 estados electrónicos, 1025 nodos fonónicos
Lobatto y cuadratura interna de orden 2. No se repitió la dinámica ni se alteró
el kernel. El diagnóstico independiente tardó **3,65 s en Geminga**, con un hilo.

| Magnitud | Mayor error relativo | Umbral original |
|---|---:|---:|
| RHS electrónico, norma L1 ponderada | 2,9703×10⁻⁴ | 10⁻³ |
| Potencias de dispersión y recombinación | 7,3653×10⁻⁴ | 10⁻³ |
| Medida fonónica sobre funciones de prueba 17/33 | 3,8811×10⁻⁴ | 10⁻³ |
| Refinamiento de las referencias | 2,5353×10⁻⁷ | 2×10⁻⁴ |

Los controles abarcan los tres campos registrados —normal, BCS con Δ=0,72 y
Γ=0, y sin gap con Δ=0,35 y Γ=0,65— y los tres perfiles térmico/no térmicos
del contrato. La mayor diferencia del RHS corresponde a dispersión del perfil
no térmico normal; la mayor diferencia de potencia corresponde a dispersión
del perfil no térmico sin gap.

La referencia integra B.7 a energía fija y luego las mismas funciones de prueba
electrónicas. Una segunda integral independiente en coordenadas (E,Ω) verifica
potencias, tasas y medidas fonónicas; sus intervalos terminan en cada nudo de
la base de 33 funciones. La base de 17 se obtiene mediante la identidad exacta
entre estas dos bases anidadas. Ninguna referencia enumera los eventos del
kernel experimental. La normalización de prueba es α²F=0,03(Ω/4)² y prefactor 1;
los errores relativos son invariantes al multiplicar uniformemente ese factor.

También se aclaró el piso de error BCS observado previamente. Para el perfil
«phonon_bubble» y el canal de dispersión, las actividades puntuales difieren
1,9475×10⁻⁴ de la referencia ideal η=0, pero sólo **2,7320×10⁻⁷** de la referencia
analítica con el mismo η=10⁻⁸. El operador completo difiere 5,0408×10⁻⁵ de esta
última. Por tanto, aquel piso se debe al regulador finito y no se presenta como
error de cuadratura que desaparece al añadir estados.

Este resultado cierra la comparación estática que faltaba al mismo corte.
**No determina por sí solo la aceptación temporal ni el dictamen global de la
etapa 2.** La referencia ideal y la referencia con η finito permanecen separadas
en `selected_continuum.json`, con todas las filas y hashes de procedencia.

Comando reproducible, desde `/home/jdiaz/pysnspd`:

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  sandbox/stage2_cells/resume_20260921/close_selected_continuum.py
```

El script sólo utiliza las medidas estáticas existentes. Se ejecutó una vez
bajo `timeout --signal=TERM --kill-after=5s 240s`; no requiere un transiente.
