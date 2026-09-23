# Auditoría independiente de los lotes 4A

Se completaron los 40 controles locales y los 6 campos espaciales prescritos. La etapa 4 sigue abierta: estos lotes no avanzan el tiempo ni imponen un circuito o un borde físico dinámico. Los archivos originales se conservan intactos.

## Integridad y coste observado

El plan y los once archivos fuente registrados coincidieron con las identidades de ambas corridas al recibirlas. La versión anterior de controls.py se conserva por sus bytes exactos en previous_delivery_exact para que los nuevos diagnósticos no rompan esa trazabilidad. Se comprobaron los 52 archivos de resultados (40 JSON locales y 6 JSON más 6 NPZ espaciales), sus tamaños y SHA-256, sin ausencias ni diferencias. El lote local duró 27,23 s y el espacial 1.878,24 s (31,30 min). Los casos finos sin diferencias finitas tardaron aproximadamente 518 s cada uno.

## Qué queda acreditado

Los controles locales dieron 36 símbolos positivos y cuatro negativos, sin resultados de signo indeterminado. Los cuatro negativos corresponden al control de gradiente alto con regularizadores 0,05 y 0,10, en vacío y con población sintética. Son diagnósticos de estados excluidos; no se intentó evolucionarlos ni se reparó su signo.

Las anclas espectrales cambian la fuerza local menos de 3,13×10⁻⁸ en términos relativos; el mayor cambio del autovalor mínimo es 3,51×10⁻⁷. No hay evidencia para seguir endureciendo aquí la integración espectral. La diferencia finita cerca de amplitud 0,001 mejora aproximadamente cuatro veces al dividir el paso por dos, como corresponde a una diferencia centrada; el error absoluto final máximo es 5,66×10⁻⁸. La identidad algebraica local de potencia KWT tiene residuo máximo 4,34×10⁻¹⁹.

## Dos problemas distintos que el total ocultaba

La norma euclídea del gradiente nodal integrado no es una norma de fuerza comparable entre mallas: un nodo representa un volumen diferente al refinar. La auditoría usa Fᵢ=Gᵢ/mᵢ y ‖F‖²=Σmᵢ|Fᵢ|², donde mᵢ es la masa de cuadratura y Gᵢ la derivada integrada de energía. El máximo de corriente de arista tampoco equivale a la densidad de corriente ni a la corriente total de una sección.

| Comparación 4×2 → 8×4, grado 4, δ=0,10 | Campo suave | Núcleo suprimido |
|---|---:|---:|
| Cambio relativo de energía | +0,0000569 % | −0,2025 % |
| Cambio de norma L² de fuerza interior | +0,0506 % | −7,983 % |
| Cambio del calor KWT interior, par Korzh | +1,499 % | −28,907 % |
| Factor del calor KWT en bordes | 2,0004 | 1,9989 |

En el campo suave, el 99,919 % del calor KWT fino viene de los nodos del borde. La energía y la fuerza interior ya son estables, pero el supuesto calor total se duplica porque se han dejado mover libremente nodos con cargas de reacción de borde. Eso es una respuesta algebraica sin restricciones, no una disipación física ya admitida. Separar los bordes en este análisis diagnostica el problema; no implementa por sí solo una condición de contorno.

El núcleo estrecho conserva además un error espacial interior apreciable, aun separando los bordes. Por tanto, no basta con corregir los bordes ni con observar que la energía global cambia poco.

## El signo negativo del centro requiere resolver primero el gradiente

El caso espacial δ=0,05 grueso tiene un solo nodo negativo: el centro (80 nm, 0), con |Δ|/Δ₀=0,04. Su autovalor mínimo es −0,784856 y la incertidumbre local 2,85×10⁻⁵. El signo es claro para ese gradiente discreto. Sin embargo, ese gradiente no representa todavía bien el campo continuo prescrito.

En coordenadas adimensionales centradas, el perfil tiene fase φ=0,07x+0,025 sin(y/5) exp[−(x/10)²]. En el centro, la derivada de su amplitud es cero; la regla del producto da ∂ₓz=i·0,04·0,07 y ∂ᵧz=i·0,04·0,005. Estas derivadas analíticas permiten comprobar el mismo cierre local sin otra simulación espacial.

| δ | Γ analítico central | Γ de malla 4×2 | Autovalor mínimo con derivadas analíticas |
|---|---:|---:|---:|
| 0,05 | 0,000425217 | 0,046679683 | 1,570796 |
| 0,10 | 0,0000531205 | 0,005831491 | 1,570796 |
| 0,20 | 0,00000413039 | 0,000453429 | 1,570796 |

Las tres anclas analíticas son positivas. En δ=0,10 la malla fina baja Γ central a 0,000140506, todavía 2,65 veces la referencia analítica. El error relativo L² ponderado del campo Γ completo baja del 22,25 % al 4,26 %. El negativo grueso de δ=0,05 es, por tanto, candidato a artefacto por falta de resolución; no prueba un núcleo continuo inestable y no justifica aumentar δ para esconderlo. Los cuatro controles locales negativos independientes mantienen su interpretación original.

## Siguiente paso útil

El posprocesamiento complementario boundary_response_review.json ya aplica restricciones estáticas explícitas a los mismos seis campos y comprueba su potencia de reacción, sin nuevas consultas espectrales. Esa corrección algebraica debe conservarse en el siguiente ensayo. Después resolver mejor el gradiente del núcleo, incluyendo δ=0,05 y δ=0,10, con la referencia analítica del perfil. Conviene comparar fuerzas y calor interiores ponderados, además de energía y Γ; no adoptar un regularizador sólo porque vuelve positivo un ensayo grueso. Con esos resultados se podrá pasar a un transiente débil controlado sin fotón y, posteriormente, al núcleo dinámico.

Una evaluación 16×8 de grado 4 contiene 2.145 nodos. Escalar únicamente el coste observado de consultas espectrales sugiere unos 33 min por caso sin diferencias finitas; es una estimación, no una corrida realizada por esta auditoría. Los comandos largos deben quedar a cargo del usuario.

## Reproducción

`audit_completed.py` comprueba los hashes registrados, reconstruye las normas y particiones desde los NPZ, compara Γ con las derivadas analíticas y evalúa las tres anclas centrales. Acepta `--raw-root` y exige un `--output-dir` nuevo; no reescribe datos anteriores. Los valores completos, sus definiciones y la identificación de las fuentes están en `independent_audit.json`.
