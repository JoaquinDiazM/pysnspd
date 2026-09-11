# Modelo pySNSPD: documentos 0.2

Revisión del 8 de septiembre de 2026. Este conjunto actualiza los cuatro documentos 0.1 y añade diez figuras originales, datos de verificación y cálculos pequeños reproducibles. No modifica los módulos de producción de pySNSPD.

## Lectura

1. **A. Microscopía y funcional común:** propósito energético, derivación de [F] hasta A.6, variación espectral, envolvente, justificación de corriente y comparación uniforme.
2. **B. Cinética, energía y temperaturas:** reacciones, balances, energía–temperatura y poblaciones sintéticas de igual energía.
3. **C. Condensado, corriente y señal:** formulación conservada, cambio de coordenadas, relajación ilustrativa y circuito con resistencia prescrita.
4. **D. Síntesis, plan y verificaciones:** qué cambia, qué se verificó y cómo continuar el desarrollo.

Los Markdown en esta carpeta son editables. Los PDF compilados están en `output/pdf/modelo_v0_2` desde la raíz del repositorio; el paquete ZIP conserva esa estructura. Las imágenes usan rutas relativas a `figuras`, de modo que acompañan a los Markdown al moverlos.

## Cambios que conviene revisar primero

- A.0 y A.5.1 reconocen la justificación de continuidad del término correctivo de Vodolazov y la memoria. La energía común se propone como una ampliación a contrastar.
- A.2.2 desarrolla la traza Nambu/espín, la resta normal y la renormalización de la ecuación (20) de Virtanen et al. hasta A.6. A.2.3 escribe explícitamente el teorema de la envolvente y sus condiciones.
- Se escribe `|Delta|` mediante notación matemática en lugar del alias de amplitud de la revisión anterior.
- B.5 y C.4.3 distinguen invertir una energía a espectro fijo de sustituir una distribución por una térmica.
- C mantiene la propuesta dinámica. Las adiciones matemáticas explican el balance y explicitan fuentes/condiciones de borde; sus ejemplos no predicen latencias del detector.
- Las figuras están incorporadas físicamente a los PDF. Las diez figuras se entregan en PNG y PDF vectorial; A y D reutilizan algunas para mantener contexto.

## Cálculos incluidos

Los scripts se encuentran en `sandbox/model_v0_2` desde la raíz. Requieren Python, NumPy, SciPy y Matplotlib. No necesitan instalar pySNSPD ni descargar catálogos de producción. Para un entorno científico ya preparado:

```sh
python sandbox/model_v0_2/checks_a.py
python sandbox/model_v0_2/checks_b.py
python sandbox/model_v0_2/checks_c.py
```

Los scripts determinan las rutas desde su propia ubicación. Por ello hay que conservar las carpetas `docs` y `sandbox` al mismo nivel. Escriben figuras en `docs/modelo_v0_2/figuras` y datos en `docs/modelo_v0_2/verificaciones`.

| Script | Problema | Entorno de la ejecución registrada |
|:--|:--|:--|
| `checks_a.py` | Usadel uniforme, tres cierres de corriente, energía y derivadas | Geminga; Python 3.10.20, NumPy 2.2.6, SciPy 1.15.3, Matplotlib 3.10.8; un hilo |
| `checks_b.py` | Reacciones, cuadraturas BCS, poblaciones sintéticas y límite normal | Windows local; Python 3.12, NumPy 2.5.3, Matplotlib 3.11.1 |
| `checks_c.py` | Inversión BCS, celda termodinámica ilustrativa, circuito prescrito | Windows local; Python 3.12, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.1 |

El cálculo A conserva colas de Matsubara hasta orden inverso cuarto y comprueba convergencia entre 400 y 3200 términos. Sus controles se limitan a los puntos declarados. B y C incluyen refinamientos y residuos en sus JSON. Los ejemplos sintéticos no se describen como resultados materiales ni como transientes del detector.

En Geminga se utilizó una carpeta aislada dentro de `/home/jdiaz/pysnspd/sandbox/model_v0_2`; no se lanzaron las etapas PRE, SS ni photon. No se cambió el entorno científico remoto existente.

## Compilación de documentos

Requiere Pandoc y Tectonic. Se verificó con Pandoc 3.9 y Tectonic 0.15.0. También puede localizarse Pandoc a través del paquete `pypandoc_binary`.

```sh
python sandbox/model_v0_2/build_documents.py --compile
```

Sin `--compile`, solo se regeneran las fuentes LaTeX en `docs/modelo_v0_2/latex`. El conversor comprueba la existencia de cada figura. Tectonic puede descargar los paquetes tipográficos que falten en su caché la primera vez; requiere acceso de red en ese caso. Las fuentes LaTeX finales se entregan también para permitir otra instalación de TeX.

La generación incluyó control visual. Para esta entrega se revisaron las 50 páginas finales y sus figuras y ecuaciones; los archivos compilaron sin avisos de desbordes, símbolos ausentes o referencias indefinidas. El registro estructural y los hashes están en `verificaciones`.

## Fuentes y trazabilidad

Los documentos 0.1 proporcionados se conservaron sin sobrescritura. Se consultó la copia `memoria_02.pdf` de Geminga (173 páginas), en particular anexo C.2, ecuaciones C.4–C.10, pp. impresas 121–123. Las referencias primarias de Vodolazov y Virtanen et al. se enlazan desde A y D, con las ecuaciones y versiones utilizadas. Los PDF de terceros no se redistribuyen en el paquete.

La revisión local del código usada como contexto es `391bb301a6c36f4438754a51a10a7fe2aade0555`. El identificador histórico citado en 0.1 se conserva únicamente como procedencia de aquella inspección. El manifiesto identifica los archivos nuevos; no atribuye al código actual todas las simulaciones históricas de la memoria.
