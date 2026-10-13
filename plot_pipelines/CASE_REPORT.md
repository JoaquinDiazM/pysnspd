# Case report individual del modelo efectivo

Pipeline: /home/jdiaz/pysnspd/plot_pipelines/07_case_report.py

Produce **un PDF autónomo por run** del formato efectivo de pipeline05
(pysnspd.effective_photon_run.v1). Acepta varios directorios como comodidad de
ejecución: nunca superpone, resta ni compara casos. Los holds sin fotón y los
formatos históricos no se presentan con el reloj de una corrida fotónica.

## Contenido estándar

1. Parámetros exclusivos del caso: geometría, espesor y resolución separados;
   material, corriente DC, absorción/deposición, energía fotónica y retenida,
   calibradores Allmaras, circuito, discriminador y precisión numérica.
2. Observables de laboratorio: Vout menos basal DC, corriente de lectura
   Ib−Is, cruce, confirmación, latencia relativa, final guardado y pendiente
   local. También Ib, Is, Vdev, vc y potencia en la carga.
3. Estado interno: curvas de gap, energías y temperaturas equivalentes;
   snapshots del gap completo y de la región central; fase relativa al DC;
   Te, Tph, exceso energético, Js,x, Jn,x, |Jtotal| y potencial.
4. Balance independiente, precisión, procedencia, hashes y límites.

Con seis snapshots automáticos son **16 páginas A4 apaisadas**. Los mapas
comparten escala de color entre los tiempos del mismo campo. Se muestran
unidades físicas: nm, ps, meV, K, kJ/m³, A/m², V y A con prefijos legibles.

Vout no incluye filtro, amplificador, ruido ni ancho de banda experimental.
Su máximo es el máximo **observado**, no el pico de un pulso completo.
La pendiente no predice jitter. El reloj empieza en la transferencia efectiva
postcascada, no en la absorción óptica absoluta. Sólo un cruce confirmado
produce latencia; un caso censurado queda sin latencia.

## Comando para cuatro casos terminados

Ejecutar en geminga. Cada repetición requiere un destino nuevo.

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/nominal775 \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/bias065_775 \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/bias075_775 \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/nominal1550 \
  --output-dir /home/jdiaz/scratch/case_reports/bias_20261013_01 \
  --workers 2
~~~

Para otro caso, sustituir la lista de --run-dir; uno solo también funciona.
--workers 2 paraleliza dos PDFs independientes, con un núcleo físico y
bibliotecas de un hilo por trabajador. Revisa afinidad, cuotas y memoria,
reserva núcleos completos y respeta el techo 90 %. Los dos ejemplos nominales
tardaron unos seis segundos por caso; no es una garantía para archivos mayores.

## Caso aún en ejecución

Añadir --allow-incomplete con un destino nuevo. Ejemplo:

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/bias075_1550 \
  --output-dir /home/jdiaz/scratch/case_reports/bias075_1550_partial_20261013_01 \
  --allow-incomplete --workers 1
~~~

Captura una copia coherente de un checkpoint atómico; no detiene ni escribe
en la simulación. Todas las páginas quedan PROVISIONALES si no hay terminal
completo. El final de trayectoria y el último snapshot se indican por separado.
No convierte ausencia de cruce hasta el checkpoint en censura a 200 ps.
Si el caso ya terminó, usa y autentica sus productos terminales.

## Elegir instantes

--snapshots-ps 0 1 5 10 15 20 solicita las muestras guardadas más cercanas;
se incluyen siempre inicial y final. Máximo 24 instantes pedidos, dentro del
intervalo de snapshots existente. No se interpola dinámica en tiempo.

Sin esa opción, una detección muestra inicial, 1 ps, primera hotbelt guardada,
muestras próximas a cruce/confirmación y final. Sin detección: inicial,
1 ps, 10 ps, cuarto, mitad y final de la ventana, eliminando duplicados.
La hotbelt es conectividad del grafo nodal con |Delta|≤0.2|DeltaDC|
entre bordes transversales; no es el gatillo eléctrico ni prueba continuidad
entre muestras.

## Datos originales y reconstrucciones

Las curvas y mapas de gap/energía/temperatura se leen de trajectory.npz y
snapshots.npz, autenticados por manifest.json, junto con el DC y la malla
original. Se rechaza corrupción o una versión distinta de las tres fuentes
utilizadas para reconstruir corrientes.

El gap se interpola como campo complejo P1 antes de tomar su magnitud.
La fase se muestra en nodos, sin interpolar ángulos y ocultando amplitudes
menores a 0.2 del DC. Un salto de fase envuelta no prueba un vórtice.

Los flujos superconductores se recuperan de la derivada gauge exacta del
gradiente energético, independiente de las ocupaciones no guardadas.
El potencial y flujos normales se recuperan de la misma neutralidad algebraica,
con Is del estado aceptado exactamente coincidente con cada snapshot.
Se contrasta Vdev y la conservación de corriente con los valores guardados.

J por triángulo es una **aproximación de visualización**: mínimos cuadrados
ponderados de proyecciones sobre caras duales conservativas, divididas por
espesor físico y longitud dual. El recibo conserva residuos de ese ajuste.
Las flechas indican dirección con longitud visual fija de 6 nm; se omiten
magnitudes menores al 0.01 % del máximo entre los snapshots seleccionados.
El color aporta la magnitud o componente firmada. No se interpretan cambios
locales de signo como inversión de corriente global.

Te y Tph son observadores de energía equivalente, no estados térmicos
impuestos. La energía electrónica guardada representa cuasipartículas y
excluye condensación/gradiente; el balance total sí incluye esos términos.
No se inventan distribuciones completas ni mapas de colisiones o calor KWT
que los snapshots no guardaron. No se evalúa el RHS cinético ni un nuevo DC,
no se evoluciona el sistema y no se aplica un ajuste experimental.

## Productos y validación

En el destino scratch: case_report_<caso>_<identidad>.pdf, un recibo JSON
por PDF, case_reports_manifest.json como inventario y .build para capturas.
Nunca reemplaza destinos existentes. No modifica datos de entrada.

Validación del 13 de octubre de 2026: 20 pruebas rápidas; ejemplos nominales
775 y 1550 nm, detección y censura respectivamente, 16 páginas cada uno.
Revisión visual de todas las páginas, escalas y etiquetas; cotejo de hashes
y reconstrucción de Vdev con error máximo aproximadamente 1.4e-17 V.
Los 13 archivos físicos del checkpoint operativo permanecen sin cambios.
Los ejemplos están en /home/jdiaz/scratch/case_report_validation_20261013_release.
