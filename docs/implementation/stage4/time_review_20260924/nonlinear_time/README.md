# Siguiente cálculo: evolución térmica no lineal hasta 1 ps

Se mantienen la misma malla de 65×65 nodos, 256 frecuencias, temperatura 0,9 K,
tiempos KWT, potencial normal y contactos del ensayo admitido. Cada evaluación
temporal resuelve el espectro no lineal actual. El operador J0 guardado se usa
para integrar eficientemente la parte rígida mediante ETD2; la corrección
incluye el resto de la fuerza no lineal, sin modificar el sistema físico.

Las funciones phi evitan invertir J0 y evitan evolucionar una coordenada
constante auxiliar. El control del paso mide la corrección en unidades de la
perturbación inicial, con una segunda trayectoria registrada para contrastar
desplazamiento, corriente, fuerza y torque exacto. Las colas diminutas no fijan
por sí solas una tolerancia relativa. Se conservan las derivadas de movilidad
y potencial que actúan sobre el residuo de la referencia.

El [plan](plan.json) y el [comando vigente](../../../../GEMINGA_COMMANDS.md)
especifican salidas, recursos y observación. Los trabajadores poseen frecuencias
independientes y reutilizan sus factorizaciones. Se usan como máximo 27 más el
coordinador, con un hilo BLAS/OMP por proceso y presupuesto compartido del 90 %.
La consola muestra paso aceptado/rechazado, avance físico y ETA.

El horizonte completo se deja al usuario: estimación provisional 5–30 minutos,
dependiente del número de pasos adaptativos. Se guardan checkpoints y hashes.
Un piloto con `--stop-after-ps 0.001` sólo acorta el horizonte de observación;
su finalización no acredita el transiente completo de 1 ps.

Se mide energía libre térmica a temperatura fija y disipación KWT/normal.
No se interpreta como energía interna conservada de electrones y fonones.
No hay fotón ni circuito efectivo nuevo. El circuito de la memoria permanece
como contrato del ensayo posterior del dispositivo.

## Piloto completado

El piloto de 0,001 ps terminó en **48,95 s**, con **4096 raíces espectrales**.
El máximo residual espectral fue 8,43×10⁻⁸, bajo la tolerancia 10⁻⁷; el defecto
instantáneo de potencia fue ≤5,62×10⁻¹³. Se conservaron los contactos.
Ambas configuraciones tomaron el mismo paso: esto comprueba que los componentes
funcionan juntos, pero no demuestra refinamiento temporal independiente.
El contraste completo hasta 1 ps sigue pendiente. Véanse el [recibo del piloto](pilot_receipt.json)
y las [nueve pruebas del integrador y comparación](etd2_tests_receipt.json).
