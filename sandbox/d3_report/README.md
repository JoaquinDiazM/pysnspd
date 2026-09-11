# Fuente del informe D3

`build_report.py` conserva el generador del informe histórico
`output/pdf/Informe_D3_proyeccion_energetica.pdf`. Fue trasladado desde
`tmp/pdfs/d3_report` al preparar pySNSPD 1.0.0; la raíz del repositorio se ajustó
a su nueva ubicación. El contenido y los cálculos del informe no se modificaron.

El programa necesita resultados y páginas renderizadas preexistentes que se
mantienen fuera de Git, bajo `tmp/pdfs/d3_report`:

- `D3_energy_projection_cache.npz`, con el diagnóstico del transitorio guardado.
- `assets/temporal/page-1.png`.
- `assets/colormaps/page-1.png` a `page-6.png`.
- `assets/profiles/page-1.png` a `page-3.png`.

Las imágenes proceden de los PDF de diagnóstico D3; el archivo
`D3_energy_projection_manifest.yaml`, guardado junto al caché, identifica los
datos externos de la ejecución original. Las dependencias del generador son
NumPy, Pillow y ReportLab. La composición original utiliza las fuentes Arial
de `C:/Windows/Fonts` y requiere esos archivos en Windows.

Con las entradas y dependencias disponibles, desde la raíz del repositorio:

```powershell
python sandbox/d3_report/build_report.py
```

El comando reemplaza el PDF histórico indicado. No forma parte de las pruebas
del paquete ni de la construcción del catálogo nuevo. La preparación de la
versión 1.0.0 conservó el PDF existente y no ejecutó este generador.
