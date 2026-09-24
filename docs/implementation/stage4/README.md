# Etapa 4: desarrollo cerrado; entrada a etapa 5 preparada

La [entrega vigente](closure_20260924/README.md) reúne el informe físico final,
la decisión de cierre y los datos compactos de dos referencias con corriente y
cinco respuestas acopladas. La etapa se cierra como desarrollo sin fotón;
**el contrato no térmico general y la promoción a producción siguen pendientes**.

El [informe final](closure_20260924/Informe_cierre_etapa_4.md)
([PDF](../../../output/pdf/implementation/Informe_cierre_etapa_4.pdf)) muestra
mapas del condensado, perfiles y relajación temporal, respuesta espacial,
inductancia, señal del circuito y la frontera actual de corriente/calor.

La admitancia compleja cambia 1,589 % al refinar y Vout sólo 0,0005825 %.
Esto admite el control reactivo. La disipación absoluta no se certifica:
el calor radial candidato supera la potencia de puerto en 31,37 % en el
refinado; falta un balance independiente con trabajo y reservorios.
El siguiente trabajo se concreta en la [entrada a etapa 5](../stage5/README.md).
No hay otra ejecución larga necesaria para esta entrega.

## Evidencia retenida

- [Campaña acoplada: plan y derivaciones](coupled_closure_20260924/README.md).
- [Euler KWT dual y control longitudinal](final_kwt_20260924/README.md):
  1 ps, 153,81 s, 1712 nodos y 256 frecuencias; margen práctico satisfecho.
- [Núcleos autoconsistentes](self_consistent_review_20260924/README.md).
- [Energía espacial](spatial_energy_20260924/README.md).
- [Separación de escalas](final_kwt_20260924/quasiclassical_assessment.md).
- [Archivo y restauración de informes intermedios](../ARCHIVOS_ARCHIVADOS.md).

La malla dual, el Euler KWT heredado de primer orden, el circuito completo,
la política de horizonte Vout más margen y v1.0.0 se conservan. La cinética
electrón-hueco no se sustituye por un ajuste instantáneo.
