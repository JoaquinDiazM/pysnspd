# Modelo pySNSPD: documentos 0.3

Revisión del 8 de septiembre de 2026. Cinco documentos en PDF y Markdown, con 17 figuras originales y verificaciones ligeras reproducibles. Esta entrega actualiza A–D y añade E, el cuaderno de aprendizaje. Los módulos del solver permanecen fuera de esta actualización.

## Lectura y respuestas

| Documento | Páginas | Contenido principal |
|:--|--:|:--|
| A. Microscopía y funcional común | 10 | Referencia normal derivada; funcional estático; λ y ajuste del gap de Simon |
| B. Cinética, energía y temperaturas | 19 | Ocupaciones FD, burbuja y Fano/escape; energía adiabática trasladada desde A |
| C. Condensado, corriente y señal | 18 | gTDGL de partida y forma estable de la memoria; conexión con la propuesta |
| D. Síntesis, plan y verificaciones | 6 | Mapa de ecuaciones trasladadas, resultados y próximos pasos |
| E. Cuaderno de aprendizaje | 14 | Cinco temas con figuras, ejemplos resueltos, ejercicios y criterio de evaluación |

Los PDF suman 67 páginas y están en `output/pdf/modelo_v0_3` desde la raíz del repositorio. Los Markdown editables están en esta carpeta. Para responder, editar `E_cuaderno_de_aprendizaje_v0_3.md` debajo de **Respuesta E01.1**, etc. Los identificadores E01–E05 deben conservarse para seguir el progreso.

Los cinco temas comienzan **ABIERTO**. Cada uno requiere 8/10 puntos y cumplir su condición esencial para pasar a **CERRADO / APROBADO**. La siguiente revisión incorporará las respuestas recibidas, su pauta y retroalimentación; un tema sin respuestas permanece sin nota. Los ejemplos resueltos son distintos de los ejercicios evaluados. No hay soluciones ocultas de estos ejercicios en los datos entregados.

## Cambios que conviene revisar primero

- **A.2.2:** el cálculo del metal normal muestra el conteo de electrones y huecos, la energía térmica, la capacidad y la entropía antes de escribir `f_n`. A mantiene el ámbito estático; B.11 y C.12 reciben los desarrollos trasladados.
- **A.5:** diferencia entre límite difusivo y acoplamiento débil; definición de λ y alcance del reescalamiento del gap de Simon, con versión y ecuaciones identificadas.
- **B.1, B.3.3 y B.8:** `f_FD` es una ocupación; `f_e^FD`, una densidad de energía libre. Absorción, retención, reparto Fano y escape se separan. La posibilidad de depósitos con varios lóbulos se menciona sin atribuirle una probabilidad calculada.
- **C.0:** empieza con la ecuación compleja KWT/Allmaras y su reformulación de la memoria. El término temporal estable se conserva y se conecta después con la fuerza propuesta. La ubicación cotejada en la memoria es anexo C, ecuación (C.3), página impresa 121 / PDF 152.
- **D.2:** tabla de correspondencia entre etiquetas 0.2 y 0.3. Los saltos de numeración de A preservan las referencias que no se trasladaron.
- **E:** modos y campos; Nambu/espín; cálculo funcional y envolvente; electrones/huecos desde semiconductores; modos, DOS y DFPT. Cada tema cierra con vínculos precisos a A–D.

## Cálculos y figuras

Todos los cálculos de esta entrega se ejecutaron en Geminga, en `/home/jdiaz/pysnspd/sandbox/model_v0_3`, usando el entorno científico existente: Python 3.10.20, NumPy 2.2.6, SciPy 1.15.3 y Matplotlib 3.10.8, un hilo. Los tres primeros scripts reproducen las pruebas ligeras de 0.2; el cuarto añade recursos para B/E.

```sh
python sandbox/model_v0_3/checks_a.py
python sandbox/model_v0_3/checks_b.py
python sandbox/model_v0_3/checks_c.py
python sandbox/model_v0_3/figures_learning.py
```

Los scripts encuentran la raíz a partir de su ubicación. Requieren NumPy, SciPy y Matplotlib; no requieren instalar pySNSPD ni cargar los catálogos materiales de producción. Las salidas están en `docs/modelo_v0_3/figuras` y `verificaciones`.

Hay 17 figuras en PNG y PDF vectorial: 2 de A, 5 de B, 4 de C y 6 de E. `checks_a.py` genera también C.4, porque la comparación uniforme pertenece ahora a C.12. Las siete figuras nuevas son B.5 y E.1–E.6. Las unidades, parámetros y residuos se registran en JSON/CSV; D resume su alcance. Los perfiles pedagógicos no son predicciones de una cascada espacial ni de latencia.

No se ejecutaron las etapas PRE, SS o photon de producción, ni un cálculo DFT/DFPT nuevo. Los ejemplos del cuaderno verifican álgebra y normalizaciones, no calibraciones materiales.

## Compilación y revisión

Se utilizó Pandoc 3.9 para convertir Markdown a LaTeX y Tectonic 0.15.0 en Geminga para generar los PDF.

```sh
python sandbox/model_v0_3/build_documents.py --compile
python sandbox/model_v0_3/verify_artifacts.py
```

Sin `--compile`, el primer comando genera únicamente `docs/modelo_v0_3/latex`. Si Pandoc no está en PATH, también puede encontrarse mediante `pypandoc_binary`. Tectonic puede necesitar descargar paquetes de su distribución tipográfica la primera vez. `verify_artifacts.py` requiere `pypdf` y comprueba etiquetas de ecuación, figuras y avisos de compilación.

Se inspeccionaron visualmente las 67 páginas finales; los documentos contienen sus 165 etiquetas de ecuación y 17 figuras, sin avisos de desbordes, símbolos ausentes o referencias indefinidas. Los registros están en `verificaciones/QA_estructura.json`, `QA_visual.json` y `compilacion/`. El cuaderno conserva los espacios de respuesta en Markdown; el PDF sirve como versión de lectura.

## Archivos y copia en Geminga

La entrega se copia bajo `/home/jdiaz/pysnspd` con la misma estructura:

- `docs/modelo_v0_3`: Markdown, README, figuras, fuentes LaTeX y verificaciones.
- `sandbox/model_v0_3`: scripts de cálculos, figuras, compilación, controles y empaquetado.
- `output/pdf/modelo_v0_3`: los cinco PDF y el ZIP completo.

El manifiesto `verificaciones/manifiesto_v0_3.json` contiene los hashes de los archivos distribuidos. El ZIP incluye todo lo anterior salvo su propia copia. Se verifica su integridad y la coincidencia de archivos locales y remotos antes de entregar.

Se conserva la versión 0.2. La copia de la memoria y los PDF de las referencias primarias se consultaron como fuentes y no se redistribuyen. Las observaciones de 0.3 se registran mediante su hash en el manifiesto. Las citas precisas están en A–E, incluyendo las diferencias de título/paginación entre preprints y versiones publicadas.

Revisión del repositorio usada como contexto: `391bb301a6c36f4438754a51a10a7fe2aade0555`. Los resultados de la memoria no se atribuyen automáticamente a esa revisión de código.
