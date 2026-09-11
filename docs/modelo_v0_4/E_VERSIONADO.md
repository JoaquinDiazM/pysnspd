# Versionado del cuaderno E

La edición del modelo sigue siendo **0.4**. Las revisiones pedagógicas no cambian por sí solas esa edición ni las versiones de A–D.

1. **Identidad:** cada clase conserva su código E01, E03, etc. El orden E.1, E.2, etc. puede cambiar. Una sustitución que cambia el propósito de una clase recibe códigos nuevos; el anterior se marca DESCARTADA y nunca se reutiliza.
2. **Revisión de clase:** cada clase activa tiene un contador independiente r1, r2, etc. Aumenta cuando cambia una explicación, ecuación, figura o actividad de esa clase. Un cambio de portada o índice no aumenta todas las clases. Se registra el motivo: corrección, ampliación, reorganización o sustitución.
3. **Revisión del ensamblado:** E-r01, E-r02, etc. identifica el PDF y el Markdown reunidos. Aumenta al entregar un conjunto modificado, manteniendo un mapa de revisiones de clase. No representa una versión nueva de la física del modelo.
4. **Aprendizaje:** ABIERTO, CERRADO / APROBADO y las respuestas son un registro distinto. Una edición del material no acredita una clase ni borra una respuesta. Si cambia una consigna ya respondida, se conserva la consigna y la respuesta anterior con su revisión; una actividad distinta recibe identificador nuevo.
5. **Referencias estables:** se conservan las etiquetas de ecuaciones existentes. Las nuevas incorporan el código de clase, por ejemplo E03.a o E10.f. Los números de orden de secciones no se usan como requisito de aprendizaje. El material previo recomendado sólo lista clases activas y tópicos, tengan o no clase exclusiva.
6. **Trazabilidad:** antes de sobrescribir una entrega se conserva su Markdown, LaTeX y PDF en `historial/E-rNN`. El JSON `E_registro_de_revisiones.json` identifica el ensamblado actual; el manifiesto de la entrega E vincula sus fuentes, figuras, comprobaciones y PDF mediante SHA-256.

## E-r02 · 9 de septiembre de 2026

| Clase | Revisión | Motivo |
|:--|:--|:--|
| E01 | r1 → r2 | Sistema, posición/desplazamiento, campos con componentes, configuración/estado, modos y ejemplos por descripción. Ocupaciones 0 y 1 explícitas. |
| E03 | r1 → r2 | Qué se varía, perfil como argumento de la derivada, analogía vectorial y derivación de las perturbaciones nulas en extremos fijos. |
| E04 | r1 → r2 | Bandas y ocupaciones; cuatro circuitos con ampliación microscópica; demostración del conteo de huecos; comparación p–n / N–S. |
| E05 | r1 → r2 | Distribución, intervalos y delta; diferencia entre conteo, probabilidad y ocupación; historia y límites de la cuantización. |
| E06 | r1 → r2 | Propósito de las ocupaciones; analogía Newton–Lagrange con cálculo; cierre del vínculo entre modo, estado y operadores. |
| E07 | r1 → r2 | Continuidad de notación y recapitulación autónoma de las operaciones, sin rehacer su propósito. |
| E08 | r1, DESCARTADA | Sustituida por dos clases con propósitos separados. |
| E09 | r1, nueva | DFT: densidad, variación, autoconsistencia, DOS normal y alcance del modelo superconductor. |
| E10 | r1, nueva | DFPT: respuesta, rigidez, interacción, función de Eliashberg y ejemplo cinético basado en Simon et al. |

E-r01 designa retrospectivamente el cuaderno 0.4 previo a esta revisión. Su contenido se conserva sin renombrar sus clases ni sus actividades. Los nuevos gráficos son modelos pedagógicos, no mediciones ni cálculos ab initio del material.
