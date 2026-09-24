# Auditoría independiente: continuación del núcleo 4A

Las dos corridas terminaron correctamente en 46,09 minutos. El plan y las once fuentes registradas se verificaron por sus hashes, conservando la versión utilizada de controls.py; los cuatro archivos de resultados también coinciden en tamaño y SHA-256. Esta auditoría no ejecutó nuevas consultas espectrales ni transientes.

## Resultado y decisión numérica

Los 561 nodos de δ=0,05 y los 2.145 de δ=0,10 tienen símbolo local positivo, con mínimo π/2. El único nodo negativo de la malla gruesa de δ=0,05 desapareció: se confirma que aquel diagnóstico dependía del gradiente mal resuelto del perfil prescrito. Esto no demuestra estabilidad de cualquier otro campo ni sustituye una referencia física de núcleo.

| δ=0,10: 561 → 2.145 nodos | Cambio |
|---|---:|
| Energía integrada | +0,00889 % |
| Norma L² de fuerza interior | +0,572 % |
| Calor interior KWT, par heredado | −1,380 % |
| Calor interior KWT, par Allmaras | −1,961 % |
| Calor interior KWT, par Korzh | −2,551 % |
| Error L² relativo de Γ frente al perfil analítico | 4,255 % → 0,181 % |

La derivada del perfil tiene error L² relativo de 0,349 % en la malla fina. Su Γ central todavía supera en 23,7 % un valor analítico muy pequeño: la diferencia absoluta es 1,26×10⁻⁵. Esa cifra puntual debe declararse, pero no justifica bloquear todo avance hasta anularla: las fuerzas y la disipación integradas ya permiten separar el error numérico dominante anterior de la sensibilidad del cierre.

**No se recomienda otra campaña dedicada solamente a refinar este mismo perfil.** Para el próximo diagnóstico se conserva δ=0,10 como candidato declarado. Esta decisión es práctica y acotada; los cambios entre mallas no son intervalos de confianza física ni cotas matemáticas del error.

## La incertidumbre que ahora importa es del cierre

A igual malla de 561 nodos, cambiar δ=0,10 por δ=0,05 modifica el calor interior en +21,24 % con el par Korzh, +19,28 % con Allmaras y +1,27 % con el par heredado. La energía cambia sólo +0,00332 % y la norma de fuerza interior +0,266 %. El resultado confirma por qué una energía integrada próxima no decide por sí sola la dinámica: KWT pondera de manera diferente las componentes de fuerza y su distribución local.

Estas cifras corresponden a una población sintética fija y a campos prescritos. No son una calibración material de δ ni una predicción de latencia. Tampoco autorizan elegir el regularizador menor por razones de estabilidad: ambos valores son positivos en sus nuevas mallas y la física del cierre aún necesita contraste.

Los bordes fijados tienen velocidad, calor y trabajo de reacción exactamente nulos en los registros. La respuesta interior coincide con la del mismo cálculo algebraico antes de restringir el borde. El balance de potencia se conserva a precisión de redondeo. Son condiciones de control de este ensayo; no representan aún los contactos o el circuito del detector.

## Siguiente paso útil para D.4.4

D.4.4 exige contrastar física del núcleo y de la disipación con una referencia espacial adecuada. Refinar indefinidamente el candidato no aporta esa referencia. El siguiente trabajo útil es una comparación física en el mismo régimen estadístico, acompañada después por un transiente débil sin fotón con balance de energía explícito. La etapa 4 sigue abierta: faltan eventos de fase, barreras o perfiles físicos contrastados, y la identificación material de la movilidad y de la relajación cinética.

Una referencia independiente viable es el vórtice térmico radial prescrito de C.7.2: amplitud d(r)=d_eq tanh(r/ℓ₀), fase azimutal y potencial vector nulo. Se resuelve el espectro espacial de Matsubara sin exigir que la amplitud sea una solución estacionaria. Este ensayo contrastaría la fuerza de núcleo del candidato; no sería una barrera de nucleación ni una simulación completa de la cinta de 80 nm.

Para ese contraste deben respetarse tres condiciones:

- Usar temperatura impuesta y **energía libre térmica** en ambos lados. Las fuerzas de los actuales estados sintéticos a ocupación fija pertenecen a otro régimen estadístico y no pueden compararse directamente como error físico.
- Imponer la rama regular θₙ(r)=cₙr+O(r³) en el centro. No añadir θ′ₙ(0)=0: eliminaría la pendiente física. La variable θₙ/r permite formular un borde regular; en un radio interior pequeño se puede emplear la relación de regularidad equivalente.
- Tratar el radio exterior como parámetro de ese problema auxiliar. Una condición de flujo uniforme local con q=1/R mejora la aproximación respecto a imponer el valor sin corriente; un único contraste de radios debe comprobar la fuerza en el núcleo interior. La energía de un vórtice contiene una cola logarítmica, por lo que exige el mismo radio o una sustracción explícita de la cola al comparar.

En unidades kBTc y ℓ₀, la fuerza espectral térmica utiliza 2d ln(T/Tc)+4π(T/Tc)Σ[d/εₙ−sinθₙ]. No se agrega otro laplaciano K₀: la dependencia espacial ya entra por la solución espectral. Esta referencia puede evaluar el cierre finito δ, pero no calibra por sí sola KWT ni su reparto dinámico de calor.

El JSON adjunto conserva cada número, definición, identidad y comparación. El script audit_followup.py reproduce el análisis de archivos en menos de unos segundos y conserva los resultados anteriores.
