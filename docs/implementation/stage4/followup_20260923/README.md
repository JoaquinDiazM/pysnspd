# Etapa 4: resolución del núcleo y ejecución paralela

Las dos corridas de continuación terminaron. La nueva referencia espacial
independiente encuentra una discrepancia grande en la fuerza del núcleo.
**El siguiente trabajo es diagnosticar el cierre, no lanzar un transiente
costoso ni repetir refinamientos ciegos. La etapa 4 sigue abierta.** Se conservan los datos,
fuentes y resultados anteriores; la [auditoría independiente](independent_audit.md)
explica las comparaciones y sus límites.

## Qué resolvieron las corridas

El caso δ=0,05 tiene 561 nodos con símbolo local positivo; el de δ=0,10,
2.145 nodos. El único signo negativo del núcleo grueso desapareció al resolver
mejor su gradiente. Aquel resultado no era una prueba suficiente de fallo
físico del perfil prescrito. Esta corrección no elimina los contraejemplos
locales de gradiente fuerte documentados en la revisión anterior.

| Comparación | Resultado |
|---|---:|
| δ=0,10, de 561 a 2.145 nodos: energía | +0,00889 % |
| Mismo refinamiento: norma de fuerza interior | +0,572 % |
| Mismo refinamiento: calor interior, par KWT Korzh | −2,551 % |
| Error relativo L² de Γ frente al perfil analítico | 4,255 % → 0,181 % |
| δ=0,10 → 0,05, ambos con 561 nodos: calor interior Korzh | +21,24 % |

La sensibilidad al cierre ya supera la diferencia entre las dos mallas más
finas para este observable. Una energía próxima no garantiza la misma
disipación. δ=0,10 se conserva en los registros como el candidato ensayado,
sin atribuirle una calibración material ni una incertidumbre física.

Los cuatro bordes fijos tienen velocidad, calor y trabajo de reacción nulos.
Se mantienen mediante cargas conjugadas explícitas, sin modificar después la
velocidad calculada. Son restricciones del control, no contactos ni condiciones
de borde admitidas para el detector. Las corrientes se comparan sumadas por
sección y las fuerzas mediante normas volumétricas ponderadas. Véanse los
gráficos de [gradiente del núcleo](figures/core_gradient.png) y
[calor y cierre](figures/heat_and_closure.png).

## Paralelismo comprobado

El [piloto de equivalencia](pilot/parallel_equivalence_receipt.json) ejecutó dos
casos de 15 nodos tanto en serie como con cuatro trabajadores. El tiempo pasó
de **19,5554 a 8,71745 s: aceleración de 2,243×**. Los resultados numéricos JSON
y los 33 arrays de cada caso fueron exactamente iguales. Es una comprobación
de equivalencia y coste del piloto, no una promesa de aceleración para mallas
grandes.

Las consultas espectrales independientes de todos los casos comparten una cola
global. Se reutilizan claves exactamente iguales, sin redondeo, interpolación
ni cambio de ecuaciones. La recolección de consultas y el ensamblado final
permanecen en el coordinador. Una consulta no programada detiene la ejecución;
no activa un cálculo serial oculto.

Geminga dispone de 16 núcleos físicos y 32 hilos. El máximo autorizado se traduce
en **27 trabajadores y un coordinador**, con afinidad en 14 núcleos completos;
quedan reservados dos núcleos físicos. Cada proceso usa un hilo BLAS. El
presupuesto también respeta las cuotas CPU y el margen de memoria visibles
de cgroup y no reserva más del 90 % de la memoria disponible efectiva.
Estos límites y las ausencias de información quedan registrados. El piloto
utilizó sólo cuatro trabajadores.

## Referencia física y siguiente decisión

El [control GL cercano a Tc](gl_core_reference.md) conecta el núcleo de un campo
real con una solución analítica. Como no tiene giro de fase, no valida δ ni un
vórtice a 0,9 K. Tampoco convierte los estados anteriores, con población
sintética fija, en estados térmicos.

La [referencia Usadel radial a temperatura impuesta](radial_reference_README.md)
ya se calculó con su [plan registrado](radial_reference_plan.json): 512 problemas
espectrales de borde y seis evaluaciones del candidato compartieron el pool de
27 trabajadores. El tiempo registrado fue 4,852 s. Los
[resultados](raw/stage4_radial_reference_20260923/summary.json) comparan fuerzas
de energía libre sobre el mismo perfil prescrito, a 0,9 K, en el núcleo
interior r≤4ℓ₀.

| Candidato | Diferencia L² relativa frente a Usadel, R=8ℓ₀ y 256 frecuencias |
|---|---:|
| δ=0,05 | 232,95 % |
| δ=0,10 | 158,09 % |
| δ=0,20 | 75,37 % |

La norma incluye el peso radial 2πr; no divide por fuerzas puntuales cercanas
a cero. Duplicar el corte de 128 a 256 frecuencias cambia la referencia sin
corrección de cola un 1,87 %; cambiar el radio exterior de 8ℓ₀ a 12ℓ₀ cambia
la referencia interior aproximadamente 0,000106 %. La estimación explícita
de cola tampoco elimina la discrepancia del candidato.

**El núcleo físico del candidato no queda admitido.** El menor desacuerdo de
δ=0,20 no autoriza seleccionarlo para ajustar el resultado. Corresponde
examinar el cierre de energía y su respuesta espacial antes de invertir en
dinámica de núcleo; no se cambia automáticamente la formulación física.
Este ensayo térmico no es una barrera, un evento de fase, una solución
autoconsistente del condensado ni la geometría de la cinta completa. Tampoco
calibra KWT o las tasas materiales absolutas. El transiente débil sin fotón
y sus balances siguen pendientes.

La [política de horizonte](horizon_policy.md) está autorizada para futuros
transientes fotónicos: observar hasta un cruce confirmado de Vout más un margen
registrado, con techo para los casos sin disparo. **Todavía no está implementada
en el integrador experimental.** Sólo cambia el horizonte de observación; el
circuito de la memoria y los acoplamientos se mantienen hasta la parada. Los
valores operacionales pendientes no se eligen automáticamente.

Producción no cambia ni se promueve el candidato. El banco estático GLL 2D y
la referencia radial son herramientas experimentales; no sustituyen el backend
Delaunay–Voronoi de producción. La
[procedencia numérica](../review_20260923/solver_scope.md) y el
[informe anterior](../review_20260923/Informe_avance_etapa_4A.md) se conservan.

La [descomposición y derivada independiente](radial_reference_analysis.md)
identifican el término que introduce el pico sin confundirlo con un error de
derivación. El [siguiente paso](next_step.md) prepara una revisión del cierre
a nivel de energía, con el [contrato de trabajo](next_step_contract.json).
