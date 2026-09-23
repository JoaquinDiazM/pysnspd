# Inicio de 4A: implementación y primer lote manual

## Qué cambió

`KWTMobility` recibe `tau_ee_Tc_ps` y `tau_ep_Tc_ps`; sus valores por defecto
siguen siendo 0,50/2,47 ps. `metadata()` registra lo utilizado. La matriz
conserva el cambio opuesto de movilidad radial y tangencial.

`delta_regularizer_bar` se propaga por los funcionales periódico, nodal, abierto
y mixto. Se eliminaron tres denominadores `.01` independientes del parámetro.
`RectangularSpatialFunctional` reutiliza energía, fuerzas, corrientes y D.36
con un dominio completamente 2D, sin llevar sus bordes a un único modo transversal.

La consulta electrónica usa un vacío analítico nuevo y raíces causales directas
en 0,0001≤|Δ|/Δ0≤1,2 y 0≤Γ/Δ0≤2; el punto normal permanece separado.
Es registro de soporte para controles, no admisión automática del núcleo.
No se reconstruye una tabla grande ni se extrapola el catálogo publicado.
Las poblaciones quedan fijas en coordenadas de conteo durante las variaciones.

Referencia material: Tc=8,65 K, Tb=0,9 K, D=0,5 cm²/s, R□=608 Ω, espesor7 nm
y ancho80 nm. De Einstein con DOS de un espín se deriva N0=9,1533e46 J⁻¹m⁻³,
σ=234962,4 S/m y ell0=4,6985 nm. Se preserva el carácter ajustado de la referencia.

## Evidencia ligera independiente

En Geminga: **177 pruebas y 79 subpruebas**, aprobadas en 7,48 s de pytest.
Incluyen comparación dimensional de movilidad, identidades de fuerzas y corrientes,
envolvente térmica, geometría completa2D, covariancia gauge y balance de disipación.
El log y la duración del proceso están en esta entrega.

El piloto local acotado terminó en **145,2 s**. Usó630 nodos de conteo; el
control cercano al normal tomó4,4 s y el rectángulo de15 nodos138,5 s.
En el primero, el error absoluto de la derivada pasó de1,71e-7 a4,28e-8 al
dividir el paso por dos. En el segundo pasó de1,75e-8 a4,46e-9; el defecto
Noether fue2,78e-17 y los tres balances KWT quedaron ≤5,56e-17.
El símbolo local fue positivo en los15 nodos del piloto, pero su malla sólo
sirve para medir coste y probar el recorrido. **No resuelve un núcleo físico.**

La campaña final evita repetir ocho evaluaciones de energía por cada caso:
las diferencias finitas espaciales de dos pasos se reservan a dos casos gruesos
con δ=0,10. Los demás comparan estados, todas las ramas de D.36 y movilidad.
El plan y el script exactos del piloto anterior a esa optimización se preservan
en `pilot_sources/`; su resultado no se atribuye a otros bytes.

## Lotes preparados

1. **Local:40 controles.** Seis estados, dos poblaciones y tresδ, más cuatro
   comparaciones de resolución espectral. Incluye el punto normal y el soporte
   pequeño, sin afirmar equivalencia con una nube fotónica. La respuesta KWT
   local es algebraica, aplicada a la derivada parcial de densidad a gradiente
   fijo; no es toda la fuerza de Euler–Lagrange espacial.
2. **Espacial:6 estados2D.** Perfil suave y perfil suprimido impuestos; tresδ
   y dos mallas (153/561 nodos). Rectángulo160×80 nm. Su longitud no es un rango
   admitido para propagación y su depresión central no es un núcleo estacionario
   ni un vórtice. No se reduce ninguna dirección a1D.

Cada fuerza se reutiliza para las parejas KWT heredada, Allmaras y Korzh.
La temperatura equivalente iguala energía electrónica a una FD; no sustituye
la población ni fuerza una termalización. El descenso espacial instantáneo
incluye los nodos libres de borde y potencial cero como control; no representa
un reservorio o circuito físico. No se usa DOS fonónica absoluta.

Se conserva un valor negativo o no resuelto de D.36 como diagnóstico. Un error
de soporte, raíz o ejecución detiene el lote y conserva los casos previos.
No hay reintentos, recortes de rigidez ni sobrescritura automática.

## Coste y salidas

Estimaciones conservadoras basadas en el piloto local, no mediciones de este
lote en Geminga: **5–20 minutos local; 45 minutos–3 horas espacial**. La variación
depende del estado espectral y de la máquina. Un proceso y un hilo por comando;
reserva orientativa1 GB RAM (memoria típica esperada menor) y <100 MB de resultados
por lote. No conviene llenar la máquina de procesos antes de comprobar el primero.

La terminal muestra barras, etapas, tiempo transcurrido y ETA por sección;
el promedio de tareas también orienta el tiempo restante del lote. No es garantía.
Se guardan `identity.json`, `progress.jsonl`, un `result.json` por caso y
`summary.json`; los estados2D añaden `fields.npz` para mapas de amplitud, fuerza,
temperatura equivalente, disipación y todas las ramas del símbolo.
No hace falta copiar el log al chat: basta avisar que terminó o que se detuvo.

Los directorios deben ser nuevos. Si ya existen, el corredor se detiene antes
de ejecutar: conservarlos y elegir otro sufijo. El modo sin `--execute` sólo
imprime el plan, sin física. La fase `pilot` es un diagnóstico ya realizado,
no un tercer comando requerido.

## Comandos para copiar

### Primero, controles locales

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/run_campaign.py \
  --plan docs/implementation/stage4/start_20260923/campaign_plan.json \
  --phase local --output-root tmp/stage4A_local_20260923 --execute
```

### Después, estados2D

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/run_campaign.py \
  --plan docs/implementation/stage4/start_20260923/campaign_plan.json \
  --phase spatial --output-root tmp/stage4A_spatial_20260923 --execute
```

Ambos lotes son independientes y pueden lanzarse en terminales distintas si
se desea; cada uno limita BLAS a un hilo y escribe su propia carpeta.
No se lanzaron desde el agente. La política de cómputo exige entrega manual
cuando se esperan más de cinco minutos.

## Estado al entregar

Etapa4 iniciada y campaña pendiente; cero trayectorias temporales nuevas.
El ancho fotónico sigue abierto. Producción y v1.0.0 no cambian. El reservorio
superconductor previo conservaδ=0,10: su acoplamiento a un núcleo distinto
necesitará parametrización consistente antes de usarse.
