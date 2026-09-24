# Energía térmica espacial: análisis independiente

El campo espectral espacial resuelve el desacuerdo de la fuerza del núcleo del cierre local anterior. La prueba admite este oráculo térmico experimental para campos prescritos; la etapa 4 continúa abierta.

Se completaron 2048 problemas espectrales en 60.347 s. Los 12 mapas finales coinciden con sus SHA-256 y las fuentes con la identidad de ejecución. Los archivos individuales por modo permanecen en Geminga; su contenido no se volvió a verificar localmente.

## Comparación con la referencia radial independiente

| Nodos | Error L2 de fuerza frente al BVP, N=256 | Cambio de fuerza N128→256 |
|---:|---:|---:|
| 33² | 2.90268% | 1.87434% |
| 65² | 0.71952% | 1.86966% |
| 129² | 0.17952% | 1.86873% |

El error de malla disminuye aproximadamente cuatro veces al dividir el paso por dos. Las comparaciones usan el mismo corte en la malla y en el BVP, de modo que esta convergencia no elimina la dependencia de corte. La norma cartesiana incluye el origen y utiliza áreas duales; se conserva por separado la norma histórica que lo excluía.

## Energía, fuerza y fase espectral

| Variación compacta | Derivada de energía | Trabajo de fuerza | Diferencia relativa |
|---|---:|---:|---:|
| amplitude | 6.65523643 | 6.65523581 | 9.34e-08 |
| phase | 0.683350737 | 0.683350892 | 2.26e-07 |

Estas variaciones mantienen la traza espectral de contacto fija y vuelven a resolver el interior. No se añadió una cola solo a la fuerza ni un término K0.

En el perfil asimétrico, la fase del propagador anómalo del modo más bajo difiere hasta 18.0035° de la fase local del condensado. El máximo ocurre en x,y=[0.0, 0.1875] ℓ0, con |Δ|/(kBTc)=0.326905 y |f₀|=0.221508. La fase espectral es una variable auxiliar estacionaria; no se está introduciendo otra variable dinámica cuántica.

El centro tiene Δ=0 y fuerza cartesiana finita [-0.5836253206354874, 1.998825399249111e-16] en unidades adimensionales. Allí f₀=[0.07046909679853938, 2.1449955391866086e-17]; no existe una fase local del gap con la cual comparar. Este punto se conserva en la nueva norma completa, 4.3076777, frente a 4.3062876 en el registro anterior que excluía el origen.

## Decisión y trabajo pendiente

Usar la energía espacial con campos espectrales auxiliares como referencia térmica experimental de la continuación. Su acción finita produce fuerza y corriente compatibles y converge hacia el BVP independiente. No se justifica ajustar otro regularizador del momento ni invertir más malla en el cierre local descartado.

El siguiente trabajo es obtener cualquier aceleración de la cola desde esta misma energía y conectar explícitamente el contrato disipativo/cinético. El cambio observado con el corte orienta ese esfuerzo. Estos resultados no validan la movilidad temporal, un transiente, el transporte no térmico, un fotón, el circuito, un vórtice autoconsistente ni el dispositivo Korzh completo.

El archivo analysis.json contiene las cifras completas, la reconstrucción diagnóstica de corriente, la comparación de mallas, los datos del origen y la procedencia verificable. El postproceso no ejecutó nuevos problemas espectrales.
