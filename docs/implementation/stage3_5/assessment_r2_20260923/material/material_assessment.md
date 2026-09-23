# NbN de Simon 2025: unidades, escala y uso permitido

Fecha: 2026-09-23. Investigación documental y cuadraturas de tablas; **cero trayectorias y cero evaluaciones del RHS físico**. La decisión recibida es cerrar la investigación 3.5 y preparar la etapa 4 sin fotón. No se admite aquí una preparación fotónica ni una DOS fonónica absoluta para predecir el hotbelt.

La revisión permite cerrar el problema con precisión: la tabla conserva valor como **forma espectral**, pero sus bytes públicos no identifican de manera suficiente la densidad de modos por volumen. El dato que falta para las ecuaciones es el producto `g_ph(E)=N_b F_b(E)`, donde `b` identifica la base de conteo —átomo o celda—. No es necesario exigir una nueva teoría microscópica completa para obtenerlo.

## Qué comprobamos de nuevo

Se consultaron los 14 árboles de commits del [repositorio original](https://github.com/qnngroup/proj-KE-solver), hasta `5b6bd747f80016da5ccd51db73c110a8ecc6abf6`. Ninguno contiene `nbn-a2f-ph_2.dat`; todos los que contienen NbN ofrecen el mismo blob `nbn-a2f-ph.dat`. Tampoco publican los inputs estructurales ni el programa que compuso sus tres columnas. `repository_history_audit.json` deja cerrado ese inventario; no corresponde seguir buscando el mismo archivo en el mismo historial.

La descarga nueva coincide con el **payload numérico** ya auditado en etapa 1: SHA256 `9aeea0948033d771deedae40da0cb4dc59fef80ac6ea41ac4d3a67b180efc610`. El hash `e94f1127…` del manifiesto anterior identifica la copia local que incluye una cabecera adicional, no estos bytes crudos; `material_results.json` histórico ya distingue ambos. Esa cabecera local no demuestra la unidad del archivo público. No se modificó ninguno.

La comparación algebraica de archivos hermanos proporciona un control adicional:

| Archivo público | Cabecera de columna 3 | Integral de columna 3 sobre el eje guardado |
|---|---|---:|
| Al-a2f-phdos.dat | estados/THz | 3.00074610 |
| nb-a2f-phdos.dat | estados/THz | 5.98360199 |
| TiN-prim-a2f-phdos.dat | estados/THz | 6.01790853 |
| nbn-a2f-ph.dat | ninguna | 0.70842006 |
| NbN derivado de etapa 1, cola recortada | unidades originales sin admitir | 0.71995014 |

Se conservó la primera fila de cada abscisa duplicada, como en etapa 1. **No hay una normalización universal a tres modos en este repositorio.** Las diferencias son compatibles con distintas bases de celda, pero el conteo por sí solo no identifica esas bases. Las tablas originales completas permanecen en `tmp/stage3_5_r2/material_original`; se versionan hashes, URLs y resultados pequeños.

## Qué fija el código y qué no

[solver.m, líneas 47–51 y 92–97](https://github.com/qnngroup/proj-KE-solver/blob/5b6bd747f80016da5ccd51db73c110a8ecc6abf6/solver.m#L47) usa eV y ps, escribe `h=2*pi*hbar` y convierte una frecuencia en THz mediante `E=h*nu`. En consecuencia:

`E=h nu=hbar omega`, `omega=2 pi nu`, y `1 THz = 4.135667697 meV` cuando THz significa ciclos por segundo. No se sustituye `omega` por `nu`.

Para el archivo **esperado pero ausente**, el código divide una columna rotulada estados/THz por `h`, obteniendo estados/eV. La conversión de una densidad cambia la ordenada porque debe preservar el número de modos: `F_E dE=F_nu dnu`. La función `A=alpha²F` utilizada en la integral de acoplamiento no recibe ese mismo jacobiano: se conserva `lambda=2 integral A(E)dE/E`.

Aplicar a la tabla pública cualquiera de esas etiquetas sigue siendo una hipótesis. Por ejemplo, interpretar su columna 3 como estados/meV produce 2.97747453 modos en el derivado; interpretarla como estados/THz produce 0.71995014. La diferencia es un factor 4.135667697. Que el primer resultado se acerque a tres **no demuestra** la unidad ni autoriza corregirlo hasta tres.

El código también fija `N=50 nm⁻³`, `N0=15 eV⁻¹ nm⁻³`, `tau_esc=10 ps` y `lambda=1.2`; su rama Debye multiplica la DOS por `Z1=1+lambda`. Son convenciones de ese programa, no una calibración transferable en bloque. La tabla del artículo usa `N=48 nm⁻³` y `tau_esc=15 ps`, en un escenario NbN de 5 nm y Tc=10 K. No debe mezclarse con el escenario Korzh de 7 nm y Tc=8.65 K. [SI25, tabla I y apéndice VII](https://arxiv.org/html/2501.13791v3#A1.SS5).

## Celda, composición y volumen: el factor que no podemos adivinar

BG19 estudia δ-NbN estequiométrico y variantes deficientes en nitrógeno; no es el archivo DFPT de SI25. Para δ-NbN, su tabla I da 21.67 Å³ por fórmula unidad calculados y 21.16 Å³ experimentales. Cada fórmula contiene dos átomos. Esos volúmenes corresponden a 46.15–47.26 fórmulas/nm³ o 92.29–94.52 átomos/nm³. BG19 muestra ramas acústicas y ópticas y discute ambos grupos en la DOS. [BG19, tablas I/III y figuras 5–6](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.99.104508/fulltext).

Así, el valor 48 publicado por SI25 se parece al conteo de fórmulas unidad, aunque su texto denomina `N` densidad de iones. Es una **inferencia de consistencia dimensional**, no una corrección demostrada al artículo. Una celda primitiva NbN con dos átomos tiene seis ramas; una DOS por átomo puede integrar tres. Ambas son válidas si se utiliza la densidad de la misma base. Mezclarlas cambia la capacidad fonónica por un factor dos.

Quantum ESPRESSO documenta `matdyn` en estados/cm⁻¹ y normalización `3*nat`. EPW ofrece salidas y ejemplos expresados en meV. Esto demuestra que la exportación necesita un jacobiano y una base explícitos; no determina qué transformación recibió el archivo NbN sin cabecera. [QE](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html), [ejemplo oficial EPW](https://docs.epw-code.org/tutorials/tutorial_04/index.html).

## La condición mínima de admisión, sin exigir microscopía exacta

Conviene expresar la entrada material mediante dos funciones:

`A(E)=alpha²F(E)` y `g_ph(E)` [modos/(energía·volumen)].

Entonces `U_ph = integral E g_ph(E) n(E) dE`. La parte fonónica del kernel depende de `N0 A(E)/g_ph(E)` —con los factores numéricos de cada canal—. Si se conserva una notación por base `b`, cambiar `F_b→c F_b` y `N_b→N_b/c` deja `g_ph`, la energía y esa combinación del kernel intactas. **No necesitamos conocer por separado dos parámetros redundantes; sí necesitamos fijar su producto físico.**

Un cambio de normalización aplicado sólo a F, manteniendo N, no es una convención: cambia el modelo. Aunque se reajuste la amplitud inicial para conservar la energía fotónica, cambian las ocupaciones, la estimulación bosónica y las tasas de absorción/creación de pares. Una `lambda` correcta no detecta este error, porque no contiene F por separado.

Quedan dos vías concretas para una admisión futura; ninguna se activa en esta entrega:

1. **Exportación trazable:** obtener F con unidad, base, número de átomos y volumen asociados; convertirla directamente a `g_ph`. Bastan esos metadatos y el procesamiento que vincula las columnas, no repetir toda la DFPT.
2. **Cierre volumétrico efectivo declarado:** conservar la forma y fijar su escala contra una capacidad calorífica fonónica independiente, con muestra y rango térmico definidos. Sería un modelo efectivo nuevo, no una recuperación de la normalización de SI25. A20 proporciona como escenario acústico `v_avg=4912 m/s` y `C_ph=B_ph T³`, con `B_ph=2*pi²*kB⁴/(5*hbar³*v_avg³)`; puede anclar un control de baja temperatura, pero no certifica los modos ópticos ni la capacidad de toda la banda. [A20, pp. 13–16](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf). No se ha elegido temperatura de ajuste ni calculado un factor de escala.

La segunda vía evita esperar indefinidamente un archivo, pero debe contrastar también el conteo de modos y declarar su incertidumbre de material. No se ajustaría para hacer pasar una prueba numérica.

## Qué tasas se pueden sostener y qué requiere la cascada

La integral de la columna 2 da `lambda=1.21566534` para el original y `1.21066677` para el derivado común de etapa 1. La cercanía al parámetro 1.2 del código es un control de consistencia del acoplamiento, no una validación de la DOS fonónica. El corte negativo previamente autorizado modifica A en su soporte común; se conserva esa transformación histórica y no se repite con una nueva regla.

Con A y un eje de energía admitidos, la tasa electrónica para **ocupaciones n especificadas** tiene escala `2*pi/hbar integral A(E)...dE`; no necesita N y F por separado. En cambio, la evolución de la ocupación fonónica, su realimentación dinámica y la conversión de energía depositada en ocupaciones requieren `g_ph`. `unit_audit.json` incluye un límite algebraico de emisión en banda normal vacía, sin fonones y por encima de toda la banda: 169.918 ps⁻¹ (original), 168.884 ps⁻¹ (derivado). **Su inversa, unos 5.9 fs, no es el tiempo de cascada ni de termalización del detector**; sólo verifica dimensiones e ilustra que son magnitudes distintas.

SI25 evoluciona una burbuja fonónica después de la conversión inicial y no proporciona un operador e–e de alta energía. A20 distingue una escala e–ph `tau0=1.87 ns`, derivada de `tau_ep(10 K)=16 ps`, de la intensidad de su integral e–e; su valor nominal adimensional es aproximadamente 12.3 y la compara con intensidades artificialmente aumentadas. Ninguno de esos números justifica identificar una relajación BGK sintética con la cascada. [A20, pp. 13, 15–16 y 24–25](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf).

Para caracterizar una cascada futura, el mínimo útil es registrar conjuntamente: una regla de reparto e–e dependiente de energía; A y `g_ph` consistentes para emisión, absorción, recombinación y creación; escape hacia el sustrato; y transporte/promedio en profundidad. Puede ser un cierre efectivo contrastable. No se exige simular cada átomo, pero una sola constante `tau_ee` o `tau_ep` no determina esos repartos ni el reloj de transferencia.

## Decisiones de esta revisión

| Elemento | Estado al cerrar 3.5 | Consecuencia |
|---|---|---|
| Datos NbN y recorte de etapa 1 | Preservados, uso de forma condicionado | No retirar la forma ni declarar unidades inexistentes |
| `E=h nu`, jacobiano de la DOS | Convención identificada en el lector original | Prohibido perder el factor 2π o aplicar el mismo jacobiano a A |
| `g_ph=N_b F_b` | Contrato suficiente identificado; valor absoluto abierto | Futuro material debe aportarlo o calibrarlo explícitamente |
| Factor a tres modos / factor dos por átomo | No adoptados | Las proximidades numéricas no son calibraciones |
| χ, ancho y reloj de cascada | No determinados por esta tabla | Respetar la decisión de caracterizar la cascada y dejar el ancho abierto |
| Etapa 4 sin fotón | Preparación permitida por decisión del usuario | No equivale a validar predicciones materiales de hotbelt o latencia fotónica |

Reproducción ligera: `python sandbox/stage3_5/assessment_r2_20260923/material_unit_audit.py --output <ruta-nueva.json>`, con los nueve archivos indicados por sus URLs/hashes en `unit_audit.json` dentro del caché. Duración medida de las cuadraturas: 0.063 s. El script se niega a sobrescribir evidencia y no inicia simulaciones.
