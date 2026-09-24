# Siguiente ejecución: resolución de momentos y proyección de potencial

Esta campaña reutiliza los núcleos térmicos terminados. No repite su relajación ni simula un fotón. Resuelve 181 consultas espectrales independientes y reutiliza cada resultado para la respuesta cinética y la proyección estática de potencial.

La malla energética base contiene 31 puntos; la fina, 50 e incluye todos los anteriores. Se refina la región térmica y el borde del gap porque los 12 puntos anteriores sobreestimaban en un 18 % la integral conocida de la susceptibilidad térmica. La comparación principal es de fuerza, corriente y flujo **ponderado por energía**, integrados con la misma cuadratura. No impone una tolerancia arbitraria nueva sobre la DOS.

| Núcleo existente | Energías calculadas | η/Δ de referencia |
|---|---:|---:|
| Radial 65×65 | 50 | 0,02 y 0,01 |
| Radial 129×129 | 50 | 0,01 |
| Asimétrico 65×65 | 31 | 0,01 |

Las siete comparaciones de proyección reutilizan también los subconjuntos de 31 puntos. η sigue siendo una resolución numérica del contorno complejo, no una tasa de relajación física. Los perfiles de perturbación siguen siendo direcciones infinitesimales de prueba.

## Comando único en Geminga

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/resolution_campaign.py \
  --plan docs/implementation/stage4/self_consistent_review_20260924/resolution_campaign/plan.json \
  --output-root /home/jdiaz/scratch/pysnspd_stage4_moment_resolution_20260924 \
  --execute
```

Reserva provisional: **5–15 minutos**, hasta 27 trabajadores y un coordinador, con un hilo numérico por proceso. El inventario y las cuotas se comprueban de nuevo al iniciar cada fase; siempre se reserva al menos el 10 % disponible y núcleos completos cuando la topología lo permite. La proyección reserva 2 GiB por trabajador; las demás fases, 1 GiB, además del coordinador. La revisión previa calcula el espacio requerido sobre el volumen de destino. Muestra progreso y ETA por fase; guarda cada consulta concluida, sus hashes y los resultados integrados. No sobrescribe directorios ni reintenta silenciosamente un fallo.

Para comprobar rutas, hashes y recursos sin ejecutar ni escribir salidas, usa el mismo comando **sin `--execute`**. La campaña completa queda para ejecución manual porque la resolución adicional y el contorno más cercano al eje real pueden superar cinco minutos. El piloto de coste separado, si figura en esta carpeta, no constituye la campaña completa.

## Qué entregará

- `spectra/`: espectros retardados sobre los mismos gaps y el mismo exterior radial.
- `kinetic/`: respuestas longitudinales y de carga sobre esos espectros congelados.
- `projection/`: siete comparaciones entre respuesta completa, modo de carga omitido y proyección de potencial.
- `moment_comparisons.json`: cambios de cuadratura, η y malla; discrepancia de proyección D, corrección omitida C y variación numérica observada U.
- `identity.json`, `executed_plan.json`, `progress.jsonl` y `summary.json`: procedencia, recursos, progreso y cierre verificable de la ejecución.

U suma cambios observados y **no es una cota rigurosa**. Las fuerzas entre mallas se comparan como densidades; las corrientes se interpolan por ancho dual antes de evaluar el mismo momento. El caso asimétrico no tiene aún un control independiente de malla energética fina. Ningún resultado cierra automáticamente la etapa 4 ni admite un nuevo estado dinámico, un transiente del detector o la promoción a producción.
