# Revisión práctica de la etapa 4

**La integración terminó; el rechazo fue posterior.** Ambas resoluciones
alcanzaron 1 ps en 1089,49 s de cómputo total. Se verificaron 71 checkpoints y
55 pasos aceptados, sin repetir ecuaciones espectrales. El estado original
`TEMPORAL_REFINEMENT_NOT_MET` y su `failure.json` se mantienen íntegros.

El criterio pasó en 95/96 comparaciones. Hasta 0,1 ps ambos planes produjeron
los mismos pasos; sólo los dos tiempos posteriores comparan discretizaciones
diferentes, con 23/24 comprobaciones satisfechas. La única discrepancia no
admitida es el torque de fase secundario inducido por la sonda de amplitud
a 1 ps: diferencia 8,30465×10⁻⁸ frente a una tolerancia 3,76819×10⁻⁸, en las
unidades y norma de la acción discreta. Equivale al 1,500 % de su pequeño valor
inicial, pero al 0,004566 % del torque inicial de la sonda angular principal.
La normalización adicional contextualiza el resultado; no cambia el veredicto.

Se cierra **este ensayo térmico como desarrollo con límites** y se continúa
hacia la malla final. No se solicita otra trayectoria cartesiana ni un nuevo
integrador para cuantificar la cola extinguida. La etapa 4 completa sigue
abierta para el acoplamiento no térmico y los controles espaciales/puertos
que lo requieran. No se acredita el detector ni se inicia la etapa 5.

- [Informe ilustrado](https://github.com/JoaquinDiazM/pysnspd/blob/c80c0f8612d8fd5a2b390938ca5ed4b9d58d3caa/output/pdf/implementation/Informe_etapa_4_revision_practica.pdf).
- [Análisis reproducido](analysis.md) y [datos y definiciones](analysis.json).
- [Fuentes publicadas y reutilización del solver](published_methods.md).
- [Ruta física siguiente](physical_route.md).
- [Decisión de alcance](decision.json).
- [Criterios de figuras y próximos cálculos](../../REPORTES_Y_CRITERIOS.md).

La acción térmica ya admite un grafo general. El [adaptador de malla](dual_mesh/README.md)
conecta las áreas y caras duales existentes; sus tres pruebas pasaron. La malla
preparada para el próximo ensayo es `dual_mesh/resampled/mesh.npz`, con 1712 nodos.
Un remuestreo del contorno eliminó los avisos y el exceso de área de la primera
preparación, conservada como evidencia. Esto comprueba geometría y transferencia
del operador; todavía no es un transiente sobre esa malla.

El [puente KWT](kwt_bridge/README.md) llama al paso local real de la memoria tras
convertir fuerza, tiempo y potencial. Su prueba pasó y las tres reducciones del
paso muestran el error de primer orden esperado. No sustituye el ensayo acoplado.
La corriente también debe venir de la nueva acción: pasar solamente una fuerza
nueva al solver GL no basta.

La preparación geométrica de una cinta de 80×160 nm es un objeto de trabajo
para esa conexión, no una simulación ni una admisión del dominio de Korzh.
Su longitud de dos anchos es provisional. La anchura fotónica y las tasas NbN
siguen abiertas. Producción y v1.0.0 se conservan.
