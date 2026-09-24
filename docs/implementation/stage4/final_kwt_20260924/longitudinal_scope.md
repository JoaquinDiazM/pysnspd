# Qué acredita el control longitudinal y qué significa cerrar la etapa 4

El cometido original D.4.4 es contrastar el núcleo, las fuerzas, la movilidad y
la disipación. Su resultado puede ser un dominio de uso junto con una decisión
de reformular; no tiene por qué ser la aceptación de toda predicción del
detector. La autorización al cerrar 3.5 fue realizar controles sin fotón. No
incluyó inventar la cascada inicial ni completar la normalización material de
la tabla fonónica por un ajuste numérico.

## Un control adicional que no exige otra campaña espectral

Se prepara el equilibrio uniforme real de la misma suma de 256 frecuencias de
Matsubara, a 0,9 K, en la malla dual de 80 por 160 nm. El carácter provisional de
esa longitud se conserva. Los contactos fijan el incremento del gap y de las
distribuciones a cero. Los bordes laterales conservan la condición natural.

La perturbación espacial es el modo propio más suave del Laplaciano discreto:
`L v = lambda M v`, donde M contiene las áreas duales y L las conductancias de
las caras. Es un campo definido sobre todos los nodos. Su normalización es
`max(|v|)=1`. En este equilibrio uniforme, las derivadas espectrales y el
transporte longitudinal son funciones de ese mismo Laplaciano. Por ello una
perturbación inicial proporcional a v permanece exactamente en ese sector de
las ecuaciones **linealizadas**. No se supone que un fotón produzca ese perfil
ni que el detector pueda reducirse a una dimensión.

Esto permite evaluar un control de ocho amplitudes en lugar de volver a
resolver Usadel en cada instante. La identidad se contrasta contra el
`SpectralTangent` existente; también se compara el RHS embebido en todos los
nodos con `LongitudinalReciprocalBridge`. No es una interpolación espectral ni
un reajuste de las poblaciones a una temperatura.

El estado tiene una pequeña supresión inicial del gap (0,2 %) y un perfil
independiente de ocupación electrónica sobre siete energías. La cuadratura es
una elección explícita de control. No es una cuadratura material convergida.
El regulador espectral eta=0,01 no aparece como baño ni como tasa de colisión.

El avance usa Euler de la tangente radial KWT. Se contrasta una vez con mitad
de paso y con la exponencial de la pequeña matriz constante, sólo como
referencia temporal. La tolerancia práctica es 2 % de la escala inicial de la
perturbación; no se divide por una cola final que se extingue. Se comprueban
además el soporte de Pauli y el balance integrado de disponibilidad.

**Disponibilidad es la energía libre cuadrática de la perturbación.** Su
descenso coincide con las pérdidas KWT y de transporte. No es la energía
interna total, ni una inyección de calor en electrones o fonones. Este control
no determina una tasa material tau_kin: sus tiempos proceden de la movilidad
heredada y del transporte, con colisiones desactivadas.

## Qué tarea pertenece a cada decisión

| Tarea | Evidencia necesaria y estado de la decisión |
|---|---|
| Contrastar el núcleo y la acción espacial común | Se apoya en los núcleos térmicos autoconsistentes, su resolución espacial/espectral y fuerzas/corrientes de la misma acción. No se identifica un núcleo físico mediante un parámetro local regularizador arbitrario. |
| Verificar movilidad y disipación | La trayectoria térmica sobre malla dual y el paso KWT real son el control principal. El control longitudinal añade intercambio no térmico débil y su pérdida de disponibilidad. |
| Identificar tau_kin y reparto material de calor | Falta normalización material de tasas y cierre dinámico compatible de energía. Los dos controles anteriores no producen estos datos; no es un problema de tolerancia. |
| Acoplar gap, poblaciones, fase, potencial y calor a amplitud finita | El compromiso posterior de etapa 4 lo exige antes de declarar implementado todo el modelo no térmico. Los controles separados verifican sus sectores, pero no prueban esa unión general. |
| Acoplar contactos y circuito durante un transiente polarizado | La infraestructura de etapa 3 está disponible. Sus verificaciones instantáneas no acreditan todavía la trayectoria acoplada; debe resolverse antes de usar esa capacidad para un dispositivo. |
| Preparación fotónica, hotbelt y retardo de Korzh | Es la predicción de etapa 5. Requiere la transferencia fotónica admitida, las tasas materiales utilizadas, los bordes y el circuito completo, además de la unión dinámica anterior. |

No corresponde desplazar un término constitutivo ausente a una nota de error
numérico. Tampoco corresponde repetir el control térmico para medir una cola
secundaria con más cifras. La decisión final debe distinguir **núcleo y
disipación investigados** de **modelo completo no térmico implementado**. Un
cierre de investigación de etapa 4 puede concluir que la integración general
requiere trabajo adicional; no puede declarar cumplido ese trabajo mediante
un cambio de nombre de etapa.
