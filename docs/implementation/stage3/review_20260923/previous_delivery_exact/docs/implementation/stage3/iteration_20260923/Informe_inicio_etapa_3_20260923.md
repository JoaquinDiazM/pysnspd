# Etapa 3: implementación inicial y lote preparado

Piloto espacial estático | 23 de septiembre de 2026

## Etapa 3 iniciada: el funcional espacial

<b>Se implementó el primer bloque 3A:</b> una energía espacial discreta que produce tanto la fuerza cartesiana como la corriente de los enlaces. El flujo regularizado modifica el espectro local; sus contribuciones se incluyen en las dos derivadas. La implementación permanece experimental.

| Comprobación realizada en Geminga | Resultado |
| --- | --- |
| Pruebas del funcional, progreso y reuso | 74 pruebas aprobadas en 4,37 s |
| Piloto registrado de ocho celdas | Aprobado en 33.9 s |
| Fuerza frente a variación de energía | Diferencia absoluta 7,08e-12; presupuesto 1,80e-4 |
| Corriente frente a variación de energía | Diferencia absoluta 3,42e-13; presupuesto 9,83e-5 |
| Fase global y referencia helicoidal | Invariancias y signos comprobados |

![Tira de 360 nm, sección de 120 nm × 7 nm. Campo prescrito de amplitud 0,9 Delta0, flujo base q ell0 = 0,1 y modulación de fase de 0,1 rad. Población preparada térmicamente y luego congelada. Las corrientes usan las escalas electrónicas de referencia del catálogo.](figures/pilot_response.png)

El perfil es una prueba estática prescrita, no una solución estacionaria del detector. Por eso la corriente puede variar entre enlaces: todavía no se resolvió el potencial que impondrá continuidad de corriente total. Los bordes físicos y el circuito se incorporarán después de este primer bloque.

Las fuerzas y la corriente se verificaron con diferencias independientes de la energía, sin termalizar de nuevo las poblaciones al perturbar los campos. La prueba distingue la corriente conjugada de la discretización de una fórmula continua simplemente muestreada.

## Validación pendiente y ejecución con progreso

![Izquierda: las diferencias observadas quedan por debajo de presupuestos fijados antes de ejecutar el piloto. Derecha: el estado suave tiene signo positivo en la malla espectral candidata y el control negativo conserva su signo inestable. Estos valores no certifican todos los estados del modelo.](figures/pilot_checks.png)

El menor valor propio del piloto es 1,5708; la mayor incertidumbre por diferencias finitas del símbolo es 1,04e-6. El control de amplitud 0,6 Delta0 y q ell0 = 1 da -0,30662 en estas unidades y se rechaza. Queda pendiente el contraste espectral 630/1260 del lote completo; la incertidumbre indicada aquí no lo sustituye.

| Lote preparado para el usuario | Alcance |
| --- | --- |
| Seis casos en 8, 16 y 32 celdas | Vacío, población térmica y no térmica; gradientes débiles de amplitud o fase |
| Piloto de ocho celdas | Se reutiliza tras comprobar fuentes, condiciones y hashes |
| Precisión espacial | Objetivo del 1 % en respuestas registradas; incertidumbre y orden observado informados |
| Costo estimado | Aproximadamente 15-20 min; margen orientativo 12-25 min, un proceso CPU y menos de 1 GB de RAM estimado |
| Salida de terminal | Barra del caso y del lote, avance completado, tiempo transcurrido y ETA aproximada |
| Archivos persistentes | Resultados, estados iniciales, recibos de integridad y progress.jsonl |

<b>El lote largo queda preparado y no fue lanzado por el agente.</b> La orden exacta está al inicio de /home/jdiaz/GEMINGA_COMMANDS.md y no usa instrucciones de screen. La estimación se ajusta con el trabajo terminado y nunca se interpreta una evaluación en curso como avance aceptado.

Al terminar, basta comunicar que concluyó; los resultados se leerán desde Geminga. Su revisión decidirá el alcance admitido para continuar con bordes y reservorios. La etapa 3 completa permanece abierta: todavía no hay evolución espacial con circuito ni un pulso de detección.
