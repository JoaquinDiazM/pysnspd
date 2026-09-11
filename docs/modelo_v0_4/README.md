# pySNSPD: edición del modelo 0.4 · cuaderno E-r02

La actualización del 9 de septiembre de 2026 revisó el cuaderno a E-r02. La edición del modelo sigue siendo **0.4**. Para la publicación 1.0.0 del 11 de septiembre se actualiza en D la descripción editorial del cuaderno vigente; sus ecuaciones y resultados se conservan. Los documentos editables están en esta carpeta; los cinco PDF y el ZIP están en `output/pdf/modelo_v0_4`, desde la raíz del repositorio. El comprobante externo `verificacion_entrega.json` registra los hashes y los equipos donde se verificó la entrega actual. La versión 1.0.0 publica esta documentación sin activar el modelo actualizado en el solver.

El candidato continuo documentado en A–D sigue siendo **experimental y restringido**. Su diagnóstico detecta sensibilidad de núcleo, pérdida del signo positivo del símbolo principal en ciertos estados y una DOS fonónica que no supera la admisión de unidades y normalización. La revisión pedagógica no modifica el solver ni incorpora ejecuciones PRE, SS o transitorios fotónicos completos.

## Lectura y versiones

- **A:** funcional estático y catálogo Usadel; numeración A.1–A.28 y equivalencias con 0.3.
- **B:** decisiones cinéticas, fuente efectiva, transporte, datos y pruebas independientes.
- **C:** energía espacial, primera variación, circuito y diagnósticos del alcance físico y matemático.
- **D:** ecuaciones continuas, condiciones de borde, admisión y secuencia de implementación experimental.
- **E:** ocho clases activas: E01, E04, E06, E07, E03, E05, E09 y E10. E02 permanece retirada; E08 se retira y se sustituye por E09 (DFT) y E10 (DFPT). Hay 24 respuestas pendientes y ninguna clase aprobada sin evidencia de aprendizaje.

`E_registro_de_revisiones.json` registra la revisión del cuaderno, la de cada clase, los retiros y el estado de aprendizaje. El contenido previo se conserva en `historial/E-r01`. Las revisiones editoriales de E no incrementan automáticamente la edición del modelo ni las revisiones de A–D. La tabla inicial del cuaderno incluye **Material previo recomendado**; cada clase recupera en su propia explicación los conceptos que necesita.

## Reproducción del cuaderno actual

El Markdown `docs/modelo_v0_4/E_cuaderno_de_aprendizaje_v0_4.md` es la **fuente editable canónica** para las siguientes observaciones y respuestas. Para generar el PDF después de editarlo, ejecutar desde la raíz del repositorio, con Pandoc o `pypandoc_binary` y Tectonic disponibles:

```bash
python sandbox/model_v0_4/learning_revision/build_learning.py --tectonic /ruta/al/ejecutable/tectonic
```

`build_learning.py` parte directamente del Markdown y construye sólo E; sin `--tectonic` genera su fuente LaTeX. No reemplaza las respuestas con los fragmentos de montaje.

Para reproducir las figuras y comprobaciones numéricas de esta entrega, usar una copia del proyecto y un entorno Python con NumPy y Matplotlib:

```bash
python sandbox/model_v0_4/learning_revision/figures_foundations.py
python sandbox/model_v0_4/learning_revision/figures_e04.py
python sandbox/model_v0_4/learning_revision/checks_modes.py
python sandbox/model_v0_4/learning_revision/figures_dft_dfpt.py
```

Los generadores verifican ejemplos resueltos y producen las figuras nuevas. Los recursos que se conservan están incluidos en `figuras`. Los gráficos de DFT/DFPT son ejemplos sintéticos identificados como tales, sin ejecutar cálculos ab initio ni resolver las actividades del estudiante.

`learning_revision/assemble_learning.py` realizó el montaje editorial de E-r02 desde los fragmentos de esa carpeta. Sólo se utiliza para reconstruir históricamente esa entrega sobre una copia: volver a ensamblar podría sustituir cambios personales y respuestas posteriores. No forma parte del procedimiento habitual de construcción desde el Markdown.

`learning_v04.py` y `learning_revision/revise_foundations.py` se conservan como programas históricos. **No se deben ejecutar para reconstruir E-r02:** el primero pertenece a la entrega anterior y puede sobrescribir imágenes; el segundo es una migración inicial, no el ensamblador actual. `build_documents.py` sigue disponible para una reconstrucción general, pero la revisión local de E usa `build_learning.py` para conservar A–D.

## Cálculos ligeros de A–D

Los programas `checks_a_v04.py`, `checks_b_v04.py`, `checks_c_v04.py` y `checks_d_v04.py`, en `sandbox/model_v0_4`, no importan el solver de producción. Su validación original se ejecutó en Geminga con Python 3.10.20, NumPy 2.2.6, SciPy 1.15.3 y Matplotlib 3.10.8. Ese registro corresponde a la entrega anterior de A–D; no implica una ejecución remota nueva de E-r02.

Para reproducir esos cálculos en el entorno original:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/model_v0_4/checks_a_v04.py
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/model_v0_4/checks_b_v04.py --material /home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/model_v0_4/checks_c_v04.py
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/model_v0_4/checks_d_v04.py
```

B requiere la entrada original indicada, o una ubicación alternativa mediante `--material`. Su SHA-256 se registra en los resultados y el manifiesto; el archivo se audita sin modificarlo y no se redistribuye. La trayectoria de colisiones utiliza un espectro normal/Debye sintético. El archivo NbN no proporciona una tasa absoluta sin normalización verificada.

Los resultados JSON/CSV están en `verificaciones` y las figuras PNG/PDF en `figuras`. El manifiesto cuenta las figuras citadas en las cinco fuentes actuales; los recursos históricos no alteran ese conteo. Las ilustraciones B_01 y B_05 se conservan de 0.3. Regenerar resultados puede cambiar sus hashes: para conservar una entrega aprobada, trabajar sobre una copia y renovar después su validación y manifiesto.

## Construcción, revisión y distribución

`verify_artifacts.py` requiere pypdf y verifica etiquetas, imágenes incrustadas, glifos de sustitución y advertencias de composición. Los registros `QA_visual_A.json` a `QA_visual_E.json` documentan la inspección de las páginas renderizadas. Si cambia un PDF, su revisión visual debe renovarse antes de empaquetarlo. Los registros de compilación están en `verificaciones/compilacion`.

`package_results.py` coteja las clases y respuestas del Markdown con `E_registro_de_revisiones.json`, cuenta las figuras citadas y comprueba los hashes de los PDF revisados. Incluye los fragmentos y programas de `learning_revision`, sin `__pycache__`. Antes de reemplazar una entrega, conserva el manifiesto anterior bajo `docs/modelo_v0_4/historial/<revisión>/entrega-<hash>` y el ZIP anterior bajo `output/pdf/modelo_v0_4/historial/<revisión>/entrega-<hash>`. Los ZIP históricos quedan fuera del nuevo paquete para evitar archivos anidados que crezcan en cada revisión.

Tras la revisión estructural y visual, ejecutar:

```bash
python sandbox/model_v0_4/package_results.py
python sandbox/model_v0_4/verify_distribution.py
```

El verificador de distribución sólo necesita la biblioteca estándar de Python. Comprueba cada archivo contra el manifiesto y contra su copia dentro del ZIP. La integridad de una entrega local no acredita por sí sola su copia a otro equipo; cualquier sincronización y verificación remotas se registran por separado.

`observaciones_v0_4_resumen.txt` conserva el resumen editorial de la solicitud que originó la entrega 0.4. Las observaciones posteriores, cambios de E y reglas de revisión se registran junto al cuaderno. Las referencias científicas y sus límites se citan dentro de cada documento.
