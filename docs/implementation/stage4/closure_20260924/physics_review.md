# Revisión física independiente del cierre de etapa 4

La campaña terminó y sostiene el desarrollo del régimen débil espacial con
corriente, dos modos de distribución, fase, neutralidad y circuito de tres
estados. La admisión cuantitativa es condicional para la respuesta dominada por
la parte reactiva. No sostiene una admisión de calor absoluto ni de transientes
no lineales con fotón. Esta distinción permite cerrar el desarrollo de etapa 4
sin certificar una capacidad que los datos aún no demuestran.

El análisis reproducible está en `sandbox/stage4_core/analyze_final_closure.py`.
Sólo lee referencias, campos complejos, resultados y manifiestos; no recalcula
espectros, ajusta parámetros, modifica tolerancias ni corrige datos.

## Qué física se observa

En la cinta provisional de 160 por 80 nm, a 0,9 K y Tc=8,65 K:

| Referencia | I_DC | Mínimo del módulo de Delta / Delta_0 | L diferencial resuelta |
|---|---:|---:|---:|
| Desfase 2 rad | 3,36812 microA | 0,998493 | 0,195174 nH |
| Desfase 8 rad, base | 13,00222 microA | 0,974714 | 0,207701 nH |
| Desfase 8 rad, refinado | mismo equilibrio | mismo equilibrio | 0,212416 nH |

Al cuadruplicar el desfase, la corriente crece 3,860 veces y la máxima
supresión del gap crece unas 16,8 veces. Es el comportamiento suave esperado
de un condensado que se debilita con el gradiente de fase: la respuesta ya no
es una relación corriente-fase estrictamente lineal. No es una hotbelt, un
núcleo resistivo ni una determinación de la corriente crítica. Ambas referencias
siguen siendo superconductoras en todo el interior.

El avance de fase guardado es exactamente 2 y 8 rad. En el segundo caso el
índice de rama vale 1 porque el recorrido supera 2pi; el salto máximo en un
enlace del recorrido es sólo 0,216 rad. El diagnóstico confirma que no se
redujo silenciosamente la referencia de 8 rad a 8-2pi. No identifica por sí
solo un vórtice. La rigidez mínima del Hessiano proyectado es positiva
(0,01409 en la base refinada), lo que acredita esas direcciones de perturbación,
no un certificado de estabilidad global en todos los grados de libertad.

La inductancia reactiva dinámica -Im(Z)/omega queda en torno a 0,207 nH para
la referencia de 8 rad. El resto de los 10 nH declarados pertenece al circuito
exterior y permanece positivo. La partición impide contar dos veces la misma
inductancia. El cambio de polarización aumenta la inductancia local, coherente
con la pérdida de rigidez del condensado.

Las frecuencias evaluadas son 9,01184, 36,04735 y 144,18940 GHz. Sus periodos
son 110,965, 27,7413 y 6,93532 ps. Los valores 17,7, 4,42 y 1,10 ps mencionados
en la preparación son 1/omega; no son periodos ni duraciones de simulación.
Las respuestas periódicas pueden reconstruirse con Re[A exp(-i omega t)] y
deben presentarse como tales, sin denominarlas transientes fotónicos.

La amplitud radial RMS por voltaje de puerto adimensional cae desde 0,40158 a
0,069679 y 0,0057966 al aumentar la frecuencia. Parte de esta caída procede de
que el voltaje impone una oscilación de fase proporcional a 1/omega, además de
la dinámica radial: no debe atribuirse toda a un tiempo de relajación extraído
del material. Las poblaciones no se reemplazaron por una temperatura.

El solapamiento disipativo radial resuelto es muy pequeño: la norma relativa
a KWT no supera 3,231e-5 entre los casos. La movilidad residual tiene autovalores
positivos. Esta observación respalda que no se detectó un doble conteo radial
apreciable en estas referencias suaves; no autoriza extrapolar una movilidad
local instantánea a una hotbelt ni omitir las poblaciones en el transiente.

## Qué precisión se obtuvo realmente

El refinamiento cambia simultáneamente modos (6 a 10), regulador, paneles y
cuadratura. Por tanto no identifica aisladamente la causa de una discrepancia.
Tomando el resultado refinado como denominador:

| Observable | Cambio relativo |
|---|---:|
| Admitancia compleja | 1,5890 % |
| Impedancia compleja | 1,5678 % |
| Vout complejo | 0,0005825 % |
| Campo radial, norma L2 con áreas | 0,6439 % |
| Campo angular, norma L2 con áreas | 0,3199 % |
| Potencial, norma L2 con áreas | 0,3129 % |
| Corrientes de enlace, norma Euclídea | 1,1303 % |
| Inductancia diferencial resuelta | 2,2195 % |
| Calor radial del candidato | 0,5992 % |
| Parte real de Y / potencia media de puerto | 45,6630 % |

Los campos nodales comparados incluyen todos los nodos y usan las áreas duales;
las corrientes son coeficientes de enlace y su comparación no se presenta como
norma física de densidad de corriente. El cambio de inductancia equivale a
0,004715 nH, sólo 0,04715 % de la inductancia total de 10 nH. Se conserva el
resultado de 2,2195 %; no se reetiqueta como cumplimiento estricto del 2 %.

El cambio absoluto de Re(Y) es 0,8194 % de |Y_refinado|, pequeño para la respuesta
compleja, pero grande para el calor que se pretende calcular precisamente con
Re(Y). Éste es el motivo físico para admitir condicionalmente el sector reactivo
y no el calor absoluto. El criterio sobre Y complejo no es una vía para ocultar
el 45,7 % de variación disipativa.

La entrada AC es una perturbación de **1 microV de la fuente Vb**, detrás de
Lb=1 microH. En el caso central refinado, Vdev alcanza sólo 4,595 pV y Vout
0,220489 nV. La casi invariancia de Vout bajo refinamiento se debe en parte a
que el circuito externo filtra fuertemente esa excitación y domina la lectura;
no prueba por sí sola que todos los campos internos estén convergidos.
`coupled_fields.npz` almacena respuesta por v_dev=e Vdev/(kBTc), no por voltio
de fuente. Para reconstruir el campo al estímulo declarado se multiplica por
el Vdev complejo obtenido de Zdev*delta I_s, además de e/(kBTc).

## Neutralidad, corriente y limitación de la base

La neutralidad RMS relativa a la amplitud de potencial está entre 0,0115 y
0,0264 %. En el refinado vale 0,01874 %. Los dos terminales difieren sólo
0,02574 % de la amplitud compleja de corriente en ese caso.

Una revisión más física que el proxy local `divergence/link_current` integra
la corriente en 1477 cortes transversales de la cinta. Cada corte suma la corriente
orientada de los enlaces que lo atraviesan, sin reconstrucción vectorial ni
interpolación. Los cortes se sitúan a mitad de cada intervalo entre posiciones
x nodales distintas. La corriente es constante entre esas posiciones en esta
representación de enlaces. La máxima desviación compleja respecto de la corriente
sumada en el contacto izquierdo es:

- 0,8631 % con polarización de 2 rad.
- 5,7816 % con 8 rad, Omega=0,2 y 6 modos.
- 4,0454 % al refinar ese caso a 10 modos.
- 7,4647 % con 8 rad y Omega=0,05.
- 3,9982 % con 8 rad y Omega=0,8.

La variación RMS alrededor de la media espacial en el refinado es 0,8114 %;
tanto la media como el RMS ponderan la longitud de cada intervalo entre cortes.
La buena coincidencia entre terminales no elimina la desviación interior.
No es correcto afirmar que la continuidad local se cumple al 2 % en toda la
cinta. La base pequeña y el regulador son fuentes posibles de error, pero este
contraste combinado no separa sus causas. El residuo de la matriz proyectada
de 1e-15 no reemplaza este diagnóstico de campo completo.

## Auditoría de potencia y factor de promedio

El factor del calor radial guardado es consistente con las convenciones del
código. Definiendo E0=kBTc y la energía de la acción

    U0 = hbar E0/(4 e^2 R_sheet),  tau_D = hbar/(2 E0),
    v = e Vdev/E0,  Omega = 2 tau_D omega,

la disipación instantánea es (U0/tau_D) sum m gamma (dr/dtau)^2.
Una amplitud pico r exp(-i omega t) tiene promedio del cuadrado de su velocidad
igual a Omega^2 |r|^2/8. Dividir por Vdev_peak^2 da exactamente

    <P_rad>/Vdev_peak^2 = Omega^2 sum m gamma |r|^2/(16 R_sheet).

Por tanto no corresponde dividir el calor por dos para forzar un acuerdo.
En el caso refinado se obtuvo P_rad/V^2=2,49818e-4 W/V^2 y
P_puerto/V^2=1,90159e-4 W/V^2: el cociente es 1,31374. En el caso base era
2,40327. Al incluir ambos terminales con sus voltajes +/-V/2, la potencia
es 1,90144e-4 W/V^2; esa corrección insignificante tampoco explica el exceso.

La diferencia no es un flujo a los reservorios calculado independientemente.
La campaña lineal no evoluciona el depósito B.41 ni las poblaciones de segundo
orden que reciben calor. El balance independiente de energía interna a
amplitud finita sigue pendiente. El defecto no se debe presentar como calor
validado mediante una tolerancia relajada, ni compensarse inventando un término
de baño. La variación interna de corriente y la sensibilidad de Re(Y) también
impiden asignar ya una causa física única a esa diferencia.

## Decisión que permiten los datos

Cerrar etapa 4 como desarrollo e investigación de núcleo y disipación, con
respuesta débil reactiva condicional. Conservar las pruebas temporales Euler
KWT ya aprobadas y la malla dual. Preparar etapa 5 empezando por la unión
energética y de carga del transiente polarizado sin fotón, antes de interpretar
hotbelt o retardo. No repetir automáticamente campañas de tolerancias para
perfeccionar colas diminutas. Tampoco promover el calor absoluto, la continuidad
local al 2 % o el modelo no térmico completo por el mero cambio de etapa.

La preparación fotónica, su anchura aún abierta, las tasas materiales NbN y
la lectura que definirá el gatillo permanecen decisiones de entrada al ensayo
de dispositivo correspondiente. No se cambió producción ni el tag v1.0.0.
