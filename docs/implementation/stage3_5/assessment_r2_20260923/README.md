# Cierre de investigación 3.5 y preparación de etapa 4

El usuario eligió **«Cerrar investigación 3.5 y preparar etapa 4 sin fotón»**.
Se cierra el estudio de fuentes, rangos condicionales y límites de uso. No se
declara completo el dominio físico ni validada la transferencia fotónica de
Korzh. La [decisión](closure_decision.json) conserva explícitamente esa distinción.

## Resultados que cambian la siguiente etapa

- **127 variables, 15 familias:** [use_domain.json](use_domain.json) enlaza cada
  entrada con su dominio matemático, restricciones acopladas, observables,
  margen propuesto, límite de extrapolación y evidencia siguiente. Los rangos
  de planificación no son intervalos físicos medidos ni 127 elecciones independientes.
- **Reparto de la cascada histórica:** la lectura independiente de dos figuras
  de Allmaras acota aproximadamente 7–9 % de energía electrónica restante y
  91–93 % fonónica a 0,187 ps. El perfil fonónico por sí solo tampoco conserva
  simultáneamente sus radios del 50 % y 90 % con una gaussiana única.
  [Extracción, supuestos y límites](cascade/cascade_assessment.md).
- **NbN:** el historial público completo no recupera el archivo con sufijo
  `_2.dat` que espera el lector original. Falta la densidad absoluta de modos
  por volumen, no una reconstrucción microscópica total. La forma de etapa 1
  se conserva; no se renormaliza hasta tres por proximidad numérica.
  [Diagnóstico y dos vías suficientes de admisión](material/material_assessment.md).
- **Movilidad del condensado:** cambiar los tiempos heredados por los de las
  fuentes afecta de manera opuesta la respuesta de amplitud y fase. No equivale
  a multiplicar todos los tiempos de una trayectoria. El código debe permitir
  seleccionar explícitamente esos tiempos antes de comparar escenarios.
  [Cálculo algebraico](planning_scales.json) y [capacidades reales](stage4_readiness.md).

La fila `taukin` corrige una ambigüedad de unidades de r1: el valor sintético
0,7 es tiempo normalizado del solver; sólo representa 0,7 ps si `t_ref=1 ps`.
No se cambia el código ni se identifica ese BGK con la cascada de alta energía.

## Qué queda preparado

La [etapa 4](../../stage4/README.md) comienza por núcleo, estabilidad y disipación
en estados controlados sin fotón. Su [contrato](../../stage4/entry_contract.json)
ordena las modificaciones y limita las afirmaciones permitidas. El catálogo
deberá cubrir las amplitudes cercanas a cero que realmente se consulten; su
punto normal aislado no cubre el hueco entre cero y 0,08 Δ0.

Los requisitos dinámicos pendientes de etapa 3 siguen vigentes para los ensayos
que los utilicen. Los controles estáticos iniciales no necesitan una cascada
fotónica ni un empalme 2D–1D. Se mantiene el circuito de la memoria; su partición
inductiva se deriva de nuevo cuando cambie material o geometría.

No se ejecutó ningún transiente nuevo, no se activaron parámetros físicos y
no se preparó un comando largo. Las cuadraturas, extracción de figuras y
comparaciones algebraicas son ligeras. Etapa 4 está **preparada, no iniciada**.

## Entrega y reproducción

El [informe final](Informe_cierre_investigacion_etapa_3_5.md) y su
[PDF](../../../../output/pdf/implementation/Informe_cierre_investigacion_etapa_3_5.pdf)
resumen los resultados con dos figuras nuevas. `delivery_manifest.json` registra
los bytes entregados; `qa.json` identifica los controles realizados.

Desde la raíz del repositorio:

```bash
python sandbox/stage3_5/assessment_r2_20260923/verify_delivery.py
```

El verificador comprueba 127 entradas y la cadena histórica mediante las copias
exactas de las entradas que cambiaron. La [libreta vigente](../../../GEMINGA_COMMANDS.md)
contiene reproducción y política de cómputo. Los resultados r1 y los cierres de
etapas anteriores permanecen congelados; sus estados de «abierto» son históricos.
