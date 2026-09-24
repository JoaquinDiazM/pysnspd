# Respuesta armónica de carga con espectro fijo

Estado: **completado**. Los 24 trabajos terminaron en 21,2 s, con 24 trabajadores
y un coordinador, una hebra por proceso y dos núcleos físicos reservados.
Pasaron 61 controles. No quedan cálculos pendientes en esta campaña.

Este ensayo pregunta cuánto cambia la respuesta al permitir que el modo de carga
almacene una perturbación, en vez de resolverlo instantáneamente. Conserva los
núcleos, espectros y probes de la campaña de momentos. No aplica un fotón ni
cambia el circuito o el solver de producción.

Se prescribe el modo de energía infinitesimal con dependencia
`exp(-i omega t)` y se resuelve

`(-L_TT - i nu M) hT = L_TL hL`,

donde `M_i=m_i Re(g_i)`, `nu=omega tD` y
`tD=hbar/(2 kB Tc)=0.44152 ps`. La masa DOS corresponde al cierre cinético
adiabático con espectro y gap fijos. En el límite normal recupera la difusión
con el coeficiente físico original. En frecuencia cero recupera la eliminación
elíptica anterior. El ensanchamiento numérico eta no se convierte en un baño.

También se calcula la proyección de un solo perfil de energía,
`hT(E)=chi(E) v`, incluyendo su almacenamiento integrado. Esta coordenada `v`
no se identifica con el potencial electrostático de la memoria.

Se conservan por separado los fasores de las dos coordenadas cartesianas de la
fuerza compleja del condensado. Los momentos usan la misma cuadratura de energía:
fuerza, corriente y flujo de energía ponderado por E. La identidad conjunta de
torque/corriente y el residual con almacenamiento quedan registrados.

Frecuencias: `nu=0,1e-4,1e-3,.01,.1,1`. Las dos menores son anclas lentas;
ninguna cifra acredita por sí sola una aproximación AC uniforme. Se publican
`hbar omega/Delta`, `hbar omega/(kB Tb)` y `hbar omega/eta`. Los valores altos son
exploratorios, especialmente cerca de rasgos espectrales estrechos. `1/omega`
es una escala angular; el período es `2 pi/omega`.

El runner comprueba los hashes de las 181 consultas ya ejecutadas y admite
24 trabajos independientes: cuatro familias por seis frecuencias. Usa un solo
presupuesto de CPU/RAM, reserva núcleos completos, fija afinidad y una hebra por
biblioteca numérica. No sobrescribe resultados ni reinicia fallos. La ejecución
de diagnóstico se limita externamente por grupo de procesos.

En la perturbación angular de la malla 65², con eta/Delta=.01, la variación
respecto de la respuesta instantánea es:

| nu | Escala 1/omega | Corriente | Torque de fase |
| --- | ---: | ---: | ---: |
| .001 | 441,5 ps | 0,0258 % | 0,0422 % |
| .01 | 44,15 ps | 0,2576 % | 0,4220 % |
| .1 | 4,415 ps | 2,552 % | 4,195 % |

En la primera fila, el perfil de un solo potencial conserva errores del 20,55 %
y 77,53 %, respectivamente. En las cuatro familias, la variación dinámica de
corriente a nu=.001 queda entre 0,0256 y 0,0261 %. La reducción espectral a un
solo perfil domina el error en estas anclas lentas; conservar hT(E) como
restricción elíptica es un candidato útil para el siguiente control débil.
Esto no acredita esa eliminación para el evento picosegundo: en la tercera
fila hbar omega/(kB Tb)=1,92 y hbar omega/eta=11,34, fuera de un argumento
adiabático uniforme. No se ha identificado todavía el potencial de memoria ni
cerrado el trabajo de un gap que se mueve.

Frecuencia cero recupera las cuatro familias estáticas anteriores con diferencia
absoluta máxima 8,89e-16. El residual conjunto de torque/corriente es menor que
1,67e-16 y el residual dinámico integrado, menor que 2,22e-16. Estos controles
confirman la implementación algebraica; no equivalen a exactitud física de una
respuesta AC a frecuencias altas.

`analysis.json` certifica por SHA-256 y tamaño los 24 mapas completos (94,4 MB)
retenidos en `/home/jdiaz/scratch/pysnspd_stage4_charge_frequency_20260924`.
El compacto local conserva mapas angulares representativos para nu=0,.001,.1;
`compact_identity.json` enlaza sus bytes con las fuentes y resultados originales.
Los fasores y sus fases se mantienen separados de una latencia de detección.

El comando original **ya fue ejecutado**. Para repetirlo explícitamente debe
elegirse una salida nueva; este ejemplo utiliza el sufijo `_repeticion`:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/charge_frequency_campaign.py --plan docs/implementation/stage4/moment_review_20260924/charge_frequency/plan.json --spectra-root /home/jdiaz/scratch/pysnspd_stage4_moment_resolution_20260924/spectra --output-root /home/jdiaz/scratch/pysnspd_stage4_charge_frequency_20260924_repeticion --execute
```

Omitir `--execute` verifica entradas y recursos. `--pilot` selecciona una sola
familia 129² a `nu=.01` para medir coste; requiere otro directorio de salida.
