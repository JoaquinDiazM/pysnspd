# Investigación 3.5: primera revisión de rangos y márgenes

Fecha: 23 de septiembre de 2026. **Investigación con escenarios provisionales;
no admisión de un dominio físico completo ni lanzamiento de transientes.**

Se revisaron las 127 entradas de las 15 familias del inventario inicial.
El [registro de rangos y restricciones](range_ledger.json) distingue los valores
elegidos como referencia, las deducciones conjuntas, las propuestas numéricas y
lo que las fuentes no identifican. Ninguna caja de parámetros recibe una
probabilidad de confianza inventada. No se cambiaron coeficientes del solver.

## Decisiones ya recibidas

1. Priorizar formación del hotbelt y latencia relativa 775/1550 nm en el hilo
   de 80 nm de Korzh, antes del jitter estadístico.
2. Medir inicialmente el tiempo desde la transferencia al modelo; mantener
   separado el retardo óptico de la cascada no resuelta.
3. Usar D=0,5 cm²/s y R□=608 Ω como referencia material provisional de Korzh,
   sin convertir parámetros ajustados en mediciones. Conductividad y DOS normal
   se derivan conjuntamente para evitar mezclar escenarios.

4. Caracterizar primero la cascada y mantener el ancho gaussiano abierto.
   La equivalencia histórica de1,4–1,9 nm no se adopta.

Las respuestas y la solicitud posterior de investigar mejor la fracción
retenida y el ancho gaussiano están en [user_decisions.json](user_decisions.json).
La familia inicialmente sugerida de retención1/3–1 y ancho5–20 nm **no fue
aprobada** y no se usa como rango físico ni como plan de barrido.

## Investigación y resultados

| Ficha | Qué resuelve |
|---|---|
| [Material y cinética](material_kinetics.md) | Resistencia, Einstein y espín; tiempos distintos; hipótesis del gap; soporte electrónico insuficiente para trasladar sin cambios el espectro NbN. |
| [Preparación fotónica](photon_preparation.md) | Retención en trabajos de jitter/latencia2019–2020, pérdidas previas/posteriores y significado del reloj de transferencia. |
| [Caracterización de la cascada](cascade_characterization.md) | Datos de transferencia por color y requisitos previos a escoger el ancho. |
| [Del disco a la gaussiana](gaussian_initial_condition.md) | Por qué aparece el disco térmico y qué anchuras pueden inferirse al cambiar de perfil y de sector físico. |
| [Geometría y resolución numérica](geometry_numerics.md) | Longitud frente a ventana temporal, márgenes de borde, resolución, coste y propuestas de presupuesto de latencia. |
| [Circuito y observables](circuit_observables.md) | CM de tres estados, partición inductiva y límites de equivalencia con la lectura experimental. |

Cada ficha tiene un registro JSON de fuentes primarias y localizadores. Las
conversiones y figuras se reproducen sin espectros nuevos ni transientes:

```bash
python sandbox/stage3_5/research_20260923/geometry_scales.py
python sandbox/stage3_5/research_20260923/profile_comparison.py
python sandbox/stage3_5/research_20260923/verify_delivery.py
```

## Qué se puede preparar y qué queda abierto

Se puede preparar un escenario con procedencia, unidades y restricciones
conjuntas. Las longitudes candidatas y la precisión de latencia se registran
como propuestas para futuros ensayos, no como resultados de éstos. El circuito
de la memoria permanece seleccionado; el inductor adicional de96 nH de K20 no
se suma ni sustituye silenciosamente a los10 nH totales de ese escenario.

Antes de una afirmación cuantitativa del detector todavía hacen falta las tasas
materiales absolutas, la correspondencia entre preparación y cascada de la
muestra elegida, la identificación de movilidad/relajación efectiva y la
lectura experimental. Continúan pendientes los requisitos dinámicos de etapa3
y la física de núcleo de etapa4. No se propone rellenar esos datos mediante
valores arbitrarios que ajusten una sola curva.

La investigación sigue abierta: la decisión recibida exige caracterizar la
cascada antes de escoger el ancho. El plan está en [cascade_plan.json](cascade_plan.json). No hay cálculo largo pendiente en esta entrega.
Las fuentes y decisiones iniciales de 3.5 permanecen intactas para conservar
el manifiesto de cierre de etapa3; esta carpeta es la revisión vigente.

El [informe de investigación](Informe_investigacion_etapa_3_5_r1.md) resume los resultados y las figuras; [versión PDF](../../../../output/pdf/implementation/Informe_investigacion_etapa_3_5_r1.pdf).

La primera extracción de la referencia de cascada está en [cascade_digitization.json](cascade_digitization.json). Los radios que contienen50% y90% no corresponden a un único ancho gaussiano en aquellos dos instantes; [figura de diagnóstico](figures/04_cascade_shape.png). Es lectura del resultado histórico, no un transiente nuevo ni una medición de Korzh.
