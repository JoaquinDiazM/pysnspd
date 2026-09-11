# Etapa 1: admitir datos y construir el catálogo

Esta implementación es posterior al tag `v1.0.0` y se importa explícitamente
desde `pysnspd.experimental`. La ruta de producción basada en la memoria no
utiliza estos módulos.

La admisión física de un material y la verificación numérica del catálogo se
registran por separado. Una referencia BCS o un espectro Debye sintético permite
comprobar algoritmos; no convierte una tabla material sin unidades o sin
normalización verificables en una entrada admisible para NbN.

Los programas de este directorio generan las evidencias de
`docs/implementation/stage1`. El informe final resume resultados y la decisión
sobre el siguiente paso. Los comandos reproducibles y sus tiempos medidos se
mantienen también en `/home/jdiaz/GEMINGA_COMMANDS.md`, con copia en
`docs/GEMINGA_COMMANDS.md`.

Las corridas automáticas usan un hilo y un límite temporal de 240 segundos cuando
su duración no está bien establecida. Un trabajo que se prevea superior a cinco
minutos queda reservado para ejecución manual del usuario. La interrupción por
tiempo no se contabiliza como verificación aprobada.

El alcance de esta etapa termina en datos y tablas locales. La admisibilidad
espacial de D.36, los flujos entre gaps diferentes, las condiciones de borde,
el circuito acoplado y los transitorios completos pertenecen a etapas posteriores.

## Resultado de esta entrega

La implementación pasa 232 pruebas en Geminga (57 nuevas), en 22.41 s.
El catálogo y sus diagnósticos tardaron 27.37 s. La entrada NbN se rechaza por
evidencia material insuficiente y defectos de los datos. La interpolación Hermite
corrige el error de borde detectado, pero el regulador todavía sesga un 8.61% la
pendiente de corriente del caso de baja energía estudiado.

**No se admite todavía la promoción a la etapa 2.** El dictamen obligatorio para
interpretar los pilotos está en
[`catalog_admission.json`](../../docs/implementation/stage1/catalog_admission.json).
La API permite evaluar tablas con fines de diagnóstico; importar o cargar un NPZ
no equivale a admitirlo para una dinámica acoplada.

## Comandos principales

Desde la raíz, con el entorno `snspd` de Geminga y un hilo por biblioteca:

```bash
python -m pytest -q
python sandbox/stage1_catalog/audit_materials.py --nbn-path /home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat --verify-upstream
python sandbox/stage1_catalog/build_catalog.py
python sandbox/stage1_catalog/check_boundary_comparison.py
python sandbox/stage1_catalog/benchmark_queries.py
python docs/implementation/stage1/review/low_energy_regulator_reference.py
python sandbox/stage1_catalog/verify_delivery.py
```

Los programas de cálculo sobrescriben sus resultados en `docs/implementation/stage1`;
el constructor permite elegir otros destinos mediante `--output` y `--catalogs`.
Después de regenerar resultados debe revisarse el dictamen, recomponerse el informe
y escribirse un manifiesto nuevo; un manifiesto antiguo detectará el cambio.
El registro de comandos de Geminga contiene las variantes con límite temporal.

`prepare_report.py` usa NumPy/Matplotlib y prepara figuras y contenido.
`build_report.py` sólo compone el PDF y requiere ReportLab y Arial o DejaVu Sans.
En Geminga las dependencias de informe se instalaron aisladamente en
`tmp/stage1_report_deps`, sin modificar el entorno de simulación.
