# Intercambio longitudinal sobre la malla dual: control aprobado

La evolución de 1 ps completó el control en **0,078 s**, sobre 1712 nodos y
1646 grados espaciales libres. Utilizó un proceso y un hilo, fijado a una CPU
tras comprobar topología, afinidad, cuotas y memoria. No requiere otra corrida
larga. Los dos tests focalizados verifican la identidad con el Hessiano
espectral existente y la unión del operador cinético completo con el modo
invariante; también comprueban refinamiento temporal y soporte de ocupaciones.

| Resultado | Valor |
|---|---:|
| Error temporal, norma de disponibilidad respecto de la perturbación inicial | 0,003288 % |
| Defecto integrado de disponibilidad, dividido por A(0) | 0,006619 % |
| Reducción de error al dividir el paso por dos | 2,00024 |
| Disponibilidad restante a 1 ps | 66,6206 % de A(0) |
| Pérdida KWT acumulada | 33,3728 % de A(0) |
| Pérdida de disponibilidad por transporte | 1,3269×10⁻⁸ de A(0) |
| Ocupación electrónica mínima / máxima | 9,60×10⁻¹¹ / 0,127826 |
| Diferencia entre RHS modal y operador sobre todos los nodos | 4,26×10⁻¹⁴ relativa |

El criterio práctico anunciado fue 2 % respecto de la escala inicial. No se
divide por una señal final casi nula. Los errores temporales de Euler quedan
muy por debajo de ese margen, sin recortar ocupaciones ni agregar un baño
artificial para compensar un balance.

![Control longitudinal y balance explícito](figures/longitudinal_exchange.png)

La [figura y su pie completo](figures/longitudinal_exchange_caption.md)
identifican qué se representa en cada panel. El recuadro del gap muestra la
diferencia de las dos resoluciones en partes por millón del gap de equilibrio,
porque en el panel principal las curvas se superponen. La curva de transporte
queda próxima a cero; su valor final está escrito explícitamente. El dato es
pequeño porque éste es un control suave concreto, no porque se haya eliminado
el transporte de las ecuaciones.

Se retienen una amplitud espacial real del gap y siete amplitudes de población
con forma energética independiente. El modo espacial más suave es un sector
invariante **exacto de las ecuaciones linealizadas alrededor del equilibrio
uniforme**, comprobado contra los operadores completos. Esto ahorra resolver
de nuevo espectros que se conocen analíticamente. No implica que un dispositivo
fotónico pueda representarse por un modo ni por una dimensión espacial.

La preparación reduce inicialmente el gap 0,2 % y añade una perturbación
electrónica débil con forma no térmica. La energía se expresa en kBTc; las
energías, los pesos y eta=0,01 están declarados en el
[recibo](receipt.json). La cuadratura es un control finito, sin afirmación de
convergencia del continuo. Eta no entra como tasa de colisión ni escape.

La **disponibilidad** es la energía libre cuadrática de la perturbación, no la
energía interna total. Su descenso coincide con las pérdidas KWT y de
transporte. Esas pérdidas no se reinyectan aquí como calor en electrones o
fonones. La dinámica no determina una tau_kin material y no incorpora
colisiones, circuito, carga o fotón. El [alcance respecto del cierre de
etapa 4](../longitudinal_scope.md) distingue estas comprobaciones de la unión
general aún necesaria para un transiente no térmico de dispositivo.

Datos: `evolution.npz` conserva el modo, las dos resoluciones de amplitudes,
las disponibilidades y pérdidas integradas. `receipt.json` registra métodos,
criterios, fuentes y recursos efectivos. La exponencial matricial de ocho
dimensiones sólo es una referencia temporal independiente, no otro integrador
del detector.
