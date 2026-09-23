# Estado vigente de la etapa 3

La [iteración de empalme y acoplamiento](coupled_20260923/README.md) contiene
los resultados actuales y el siguiente cálculo. La campaña abierta terminó
sus seis casos en 271,28 s y conserva la rama uniforme. El piloto mixto 2D-1D
pasó en 78,34 s: energía común, fuerza, corriente, potencial, KWT, depósito de
calor y circuito de la memoria tienen balances instantáneos consistentes.

El acelerador espectral acepta sólo una semilla seguida del cálculo causal
original: 1,904 veces más rápido por kernel, sin cambiar energías ni fuerzas
por interpolaciones. Sus límites de campo son explícitos.

Los extremos libres dominan el calor de la hélice. Por ese diagnóstico, el
siguiente lote incluye carga prescrita de reservorio, amplitud terminal fija
y trabajo del borde. La comparación de mallas busca la respuesta con esa
carga; no pretende converger el calor del perfil incompatible con extremos libres.
El comando manual con barras y ETA está en `/home/jdiaz/GEMINGA_COMMANDS.md`.

La etapa 3 sigue abierta. Aún faltan la evolución temporal, la rama de borde
dependiente de la corriente y el transporte cinético a igual energía en la
unión. El empalme de campo admite la traza transversal uniforme. No se anuncia
una señal de detección ni una admisión de producción.

Las iteraciones [de bordes](ports_20260923/README.md),
[nodal](nodal_20260923/README.md) e [inicial](iteration_20260923/README.md) quedan
preservadas. El verificador vigente es
`sandbox/stage3_spatial/coupled_20260923/verify_delivery.py` y comprueba las
entregas anteriores usando instantáneas exactas de las entradas reemplazadas.

`README.md` y `entry_contract.json` de esta carpeta mantienen la propuesta del
22 de septiembre incluida en el cierre de etapa 2; su estado «no iniciada» es
histórico. El tag `v1.0.0` y los módulos de producción no cambian.
