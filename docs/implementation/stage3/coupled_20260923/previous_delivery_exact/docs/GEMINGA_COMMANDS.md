# Geminga: comandos vigentes

Actualizado el 23 de septiembre de 2026. Cuenta `jdiaz`, sin administrador.

## Pendiente: etapa 3B, equilibrios con extremos abiertos

**Preparado para ejecución manual; no lanzado por el agente.**

La campaña nodal terminó: 18 casos, diez estimaciones relativas aceptadas y una
energía de exceso sin certificado relativo por su pequeña diferencia. 3A puede
cerrarse para desarrollo en ese alcance. No repetir el lote nodal ni etapa 2.

Ya pasaron 168 pruebas y 23 subpruebas de los módulos abiertos y eléctricos.
El piloto abierto de 17 nodos pasó en 20,9 s. El circuito completo de la memoria
se verificó con una resistencia prescrita; eso todavía no es una señal SNSPD.

El siguiente lote comprueba longitudes 360/720/1080 nm con dos resoluciones por
longitud. Incluye el caso del piloto como nueva ejecución independiente. Registra
soporte, estabilidad, corriente y gradiente terminal, y partición fija de inductancia.
No integra un transiente ni una interfaz 2D-1D.

Tiempo previsto: **6-10 minutos**; puede variar con la máquina y las iteraciones.
La estimación por nodos a partir del piloto da unos 5,8-6,3 minutos; el margen
cubre trabajo del solver. Recursos: un proceso CPU, un hilo, menos de 1 GB RAM,
menos de 20 MB de salida esperada. No se ejecuta automáticamente por superar 5 min.

Ejecutar directamente:

```bash
cd /home/jdiaz/pysnspd
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage3_spatial/ports_20260923/run_open_batch.py \
  --registration docs/implementation/stage3/ports_20260923/open_registration.json \
  --output-root tmp/stage3b_open_full_20260923 \
  --execute
```

La barra mide casos completos y la ETA usa el coste por elemento de los casos
terminados. Antes del primer caso la ETA es desconocida; se imprimen consultas,
nodos y tiempo transcurrido. La velocidad puede cambiar entre mallas. Un caso
sólo pasa tras medir los residuos, no por la bandera de éxito del optimizador.

Salidas en `tmp/stage3b_open_full_20260923/`:

- `manifest.json`: criterios, casos, argumentos y SHA256 de fuentes.
- `progress.jsonl`: barra, tiempos y estimaciones de avance.
- Cada carpeta de caso: `reservoir.json`, `iterations.json`, `result.json`.
- `summary.json`: los seis casos y límites del alcance, si termina el lote.
- `failure.json`: excepción y contexto si se detiene.

El directorio debe ser nuevo. Ante fallo no se sobrescribe, relanza ni cambia
el algoritmo automáticamente. Indicar «terminó» o «reinicia» basta: el agente
leerá los archivos directamente, sin pedir que se copie toda la terminal.

## Comandos ligeros opcionales

Verificar la entrega y describir los casos sin ejecutar física:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_spatial/ports_20260923/verify_delivery.py
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_spatial/ports_20260923/run_open_batch.py
```

Últimas líneas guardadas del lote:

```bash
tail -n 8 /home/jdiaz/pysnspd/tmp/stage3b_open_full_20260923/progress.jsonl
```

Reproducir el control circuital independiente (segundos; elegir carpeta nueva):

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_spatial/ports_20260923/diagnose_electrical.py \
  --lk-ext-h 7e-9 --output-dir tmp/stage3_electrical_user_01
```

Los 7 nH son una entrada explícita de ese control aislado. La partición del
piloto físico es 0,166621 nH resueltos y 9,833379 nH exteriores para 360 nm y
8,63351 microA; no mezclar ambas pruebas ni atribuir el pulso prescrito al detector.

## Historial preservado y siguiente trabajo

Libreta anterior exacta:
`docs/implementation/stage3/ports_20260923/GEMINGA_COMMANDS_before_ports.md`.
Esa libreta conserva las referencias a archivos históricos anteriores.
Campaña nodal íntegra: `docs/implementation/stage3/ports_20260923/raw/nodal_campaign/`.
Informe actual: `output/pdf/implementation/Informe_avance_bordes_etapa_3_20260923.pdf`.

Después de revisar este lote siguen el empalme 2D-1D y el ensayo débil acoplado
en el tiempo con bordes, potencial, cinética y circuito de la memoria. La etapa 3
completa sigue abierta. Producción y el tag v1.0.0 permanecen sin cambios.

Todo cálculo previsto de más de cinco minutos se deja al usuario. No se usan
trabajos de fondo ni sondeos automáticos para eludir esa entrega. Las pruebas
breves inciertas tienen límite de 240 s; un timeout queda incompleto y se entrega
el comando sin restricción al usuario, sin fraccionarlo para eludir ese límite.
