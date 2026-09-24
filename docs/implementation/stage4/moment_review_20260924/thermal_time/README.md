# Siguiente ejecución: evolución térmica débil con todos los nodos

Estado: **preparada; el transiente físico completo todavía no se ha ejecutado**.
El preflight del operador real pasó. Seis pruebas unitarias y un fixture
sintético 5×5 de dos frecuencias pasaron; este último terminó en 0,51 s y
verificó el pool residente, checkpoints, evolución afín y contraste refinado.
Los recibos finales son `unit_tests_receipt_v2.json`, `dry_run_v2.json` y
`smoke_receipt_v3.json`. Se conservan los intentos de preparación anteriores:
uno comprobó la barrera de disco y otro detectó un Booleano NumPy no serializable
en el progreso, corregido antes de congelar las fuentes actuales.

Reutiliza el operador térmico ya validado y sus 256 matrices espectrales. Cada
trabajador conserva una fracción fija de las frecuencias y factoriza sus matrices
una sola vez. Las direcciones que solicita Arnoldi conservan todos los grados
de libertad espaciales; las dos perturbaciones iniciales no restringen la malla
a dos modos físicos.

Se integra `xdot=F0+J0 x`, con tiempo físico, temperatura y constantes KWT
heredadas. La referencia no se fuerza a ser estacionaria: se conserva `F0` y
la derivada de movilidad que actúa sobre su gradiente residual. Se guardan tanto
la evolución de esa base como la diferencia causada por cada perturbación.
Los bordes de gap y potencial permanecen fijados como en el ensayo de núcleo.
Este dominio cerrado no crea puertos ni reemplaza el circuito de la memoria.

La ventana de diagnóstico llega a 1 ps. La amplitud inicial es 0,001 veces las
direcciones registradas. Si la base o un estado perturbado sale de una vecindad
del 2 % del gap de referencia, la ejecución termina y guarda el estado: no lo
recorta ni cambia constantes. Ese margen declara el alcance lineal, no un rango
experimental ni una precisión del detector.

El criterio temporal mide el defecto espacial completo respecto de la respuesta,
excluyendo el fondo del gap y la coordenada constante de la formulación afín.
El defecto integrado objetivo es 0,1 %; no constituye por sí solo una cota del
error propagado. Un segundo pase usa más vectores de Arnoldi, tolerancia de
defecto reducida a la mitad e intervalos partidos. Se comparan desplazamiento,
fuerza, corriente y torque en los mismos tiempos físicos. Si la diferencia
supera el margen registrado de 0,5 % de respuesta/perturbación inicial, se
detiene la interpretación, sin iniciar campañas adicionales automáticamente.

Se publican energía libre cuadrática, disipación KWT/Joule aproximada y defecto
de truncación constitutivo. No se afirma conservación exacta ni monotonicidad
no lineal a partir de la ley afín. No se incluye un fotón, cinética no térmica
con gap móvil, evento de detección o circuito artificial.

Comando completo para Geminga, una vez recibida la entrega:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/thermal_time_campaign.py --plan docs/implementation/stage4/moment_review_20260924/thermal_time/plan.json --operator-root /home/jdiaz/scratch/stage4_thermal_weak_20260924 --output-root /home/jdiaz/scratch/stage4_thermal_time_20260924 --execute
```

Quitar `--execute` verifica los hashes, la aceptación del operador y los recursos
sin crear resultados ni lanzar el pool. La estimación conservadora es de
5–45 minutos, todavía sin medir sobre el transiente completo; depende de la
rigidez y los reinicios de Krylov. El programa vuelve a comprobar CPU, afinidad,
cuotas y RAM: como máximo 27 trabajadores y un coordinador, una hebra numérica
por proceso y dos núcleos físicos libres en el inventario actual. Imprime avance
temporal y ETA, y reserva 4 GiB de disco para resultados/checkpoints.

Una falla conserva cada estado aceptado, el registro del defecto y una causa
explícita. El directorio debe ser nuevo; no hay sobrescritura ni reanudación
implícita. El horizonte del ensayo no modifica ninguna escala física del sistema.
