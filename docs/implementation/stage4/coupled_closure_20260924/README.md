# Etapa 4: ensayo de unión dinámica con corriente y circuito

Se preparó un único comando para obtener los datos que faltan antes de decidir
el cierre. **El lote largo todavía no se ha ejecutado.** Primero construye dos
referencias estacionarias con corriente en la malla dual admitida; después
resuelve cinco casos de respuesta débil con espectro, poblaciones, fase,
potencial y circuito. No repite las campañas temporales ya aprobadas.

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/run_final_coupling.py \
  --plan docs/implementation/stage4/coupled_closure_20260924/final_coupling_plan.json \
  --output-root /home/jdiaz/scratch/stage4_final_coupling_20260924 --execute
```

El mismo comando está en `/home/jdiaz/GEMINGA_COMMANDS.md`. No requiere crear
sesiones. Muestra barras y ETA, conserva los registros y no reutiliza un destino
existente. Dos casos pueden avanzar simultáneamente; sus consultas espectrales
se reparten entre trabajadores dentro de **un presupuesto global de 28 hilos**
en el equipo actual de 32. BLAS/OpenMP usa un hilo por proceso y quedan dos
núcleos físicos libres. La reserva de RAM es conservadora y se comprueba al
iniciar cada etapa; no es un límite impuesto por el sistema operativo.

## Qué se va a calcular

| Parte | Cálculo | Para qué sirve |
|---|---|---|
| Referencias | Desfase entre contactos de 2 y 8 rad; 256 frecuencias Matsubara | Obtener campos estacionarios con corriente, sin imponer un gap uniforme en el interior. La corriente resultante se informa en amperios. |
| Polarización | Respuesta a la misma frecuencia en ambas referencias | Medir la mezcla entre amplitud, fase y desequilibrio electrón–hueco cuando circula corriente. |
| Rapidez | Tres frecuencias en la referencia de 8 rad | Comparar respuestas que varían en escalas aproximadas de 17,7, 4,42 y 1,10 ps, definidas como 1/ω. No son duraciones de transientes. |
| Refinamiento independiente | Más modos espaciales, cuadratura más fina, menor regulador y mayor corte energético | Comprobar que los observables del caso central no dependan de la resolución elegida. Si discrepan, la primera comparación no identifica por sí sola cuál de los cuatro cambios importa. |
| Circuito | Red completa de tres estados de la memoria | Calcular corriente del dispositivo y Vout sin acortar inductancias, capacidad ni resistencias. |

Los espectros y ambas distribuciones se resuelven en **todos los 1712 nodos**.
La respuesta del gap y del potencial usa una base espacial declarada de 6
modos, ampliada a 10 en el contraste. Es una aproximación de Galerkin sobre
la malla 2D, no una reducción física a una cinta 1D. Se guardan los residuos
fuera de esa base para que resolver exactamente la matriz pequeña no oculte
un error espacial. Los bordes laterales son aislantes y los contactos
superconductores cambian su fase de manera coherente con su tensión.

La inductancia exterior es la total de referencia menos la inductancia
diferencial de la región espacial resuelta, evaluada sobre la misma rama
estacionaria. No se resta de nuevo en cada instante. El punto DC asociado
al control se expresa mediante la corriente obtenida y la tensión de fuente
necesaria para sostenerla; los valores de los componentes circuitales no se
ajustan para acelerar el cálculo. Estas dos polarizaciones de control no se
presentan como corrientes medidas de Korzh.

## Evidencia ejecutada antes de entregar el comando

- Referencia uniforme con colas infinitas: **720 consultas, sin fallos de ejecución, en 5,60 s**.
  Se compararon tres reguladores, dos órdenes de cuadratura y dos cortes.
  La respuesta acoplada cambia como máximo 1,02×10⁻⁷ en norma relativa al
  comparar cuadraturas; el ensayo previo truncado llegaba a discrepar 55 %.
  La variación del regulador se informa aparte: el núcleo alcanza 2,16 % en
  algunos puntos, mientras admitancia y respuesta quedan por debajo de 1,71 %.
  No se interpreta esa pequeña sensibilidad como un fallo catastrófico.
  Las diferencias físicas y numéricas se separan en
  [el análisis de resultados](uniform_analysis.md).
- Una consulta espectral/cinética completa sobre los **1712 nodos**, con 13
  direcciones de perturbación, tardó **1,92 s**. Sus residuos espectral y
  cinético máximos fueron 4,95×10⁻¹² y 1,98×10⁻¹⁵. Es una medición de coste
  sobre fondo uniforme, no una estimación garantizada del fondo polarizado.
- El encadenamiento completo también se probó en una malla pequeña. La
  cuadratura deliberadamente gruesa daba una admitancia incorrecta; al usar
  la resolución prevista, el cálculo terminó en **28,78 s** y recuperó la
  referencia uniforme: Im(Y)=0,0226157 S. El defecto de continuidad medido
  frente a la escala de corriente de enlace bajó de 3,1 % a 0,000794 %.
  No se corrigieron los datos recortando estados ni forzando el balance.
- Las pruebas analíticas cubren difusión normal, cambios de calibre,
  normalización espectral, contactos móviles, conservación de carga,
  proyecciones modales y el circuito completo. La validación exacta de las
  ecuaciones auxiliares no se confunde con una precisión física exigida al
  detector.
- El flujo completo de la entrega, incluyendo referencia, recibo de aceptación,
  dependencias, respuesta y circuito, se ejecutó en una malla pequeña con corriente:
  **31,92 s**, 8 barridos estacionarios y 1112 energías. Su propósito fue comprobar
  la integración del programa; los dos modos de ese piloto no son una aceptación
  espacial del dispositivo. La [evidencia](pipeline_smoke/integration_receipt.json)
  conserva comando, fuentes y resultados.

Las colas de alta energía se integran con el cambio exacto `x=C/E`, con 16
nodos por cola en el caso principal y 24 en el refinado. C es la transición
entre cuadraturas, no un descarte de estados. El trabajador espacial también
se comprobó hasta E/(kBTc)=12000: conserva la caída esperada de los momentos
y no mostró ruido catastrófico tras aplicar el jacobiano del cambio de variable.

El piloto sugiere un lote del orden de **una a varias horas**; cerca de los
bordes espectrales y con corriente puede tardar más. La ETA se actualiza con
las consultas efectivamente completadas. Los cálculos largos quedan para
ejecución manual conforme a la política acordada.

## Qué permitirá decidir el resultado

Se conserva el margen práctico del **2 % en observables con señal apreciable**.
Una componente disipativa que tiende a cero no se convierte en un fallo grande
por dividir entre ese mismo cero. Se informan por separado la admitancia
compleja, Vout, residuos espaciales, dependencia del regulador y mezcla entre
las distribuciones. No se relajan una corriente con signo incorrecto, una
inductancia exterior negativa o una rama espectral no causal.

El núcleo radial microscópico se mide después de resolver la respuesta de
las poblaciones y del potencial. Se compara con la relajación KWT heredada.
El candidato suma esa relajación radial una vez; **no se adopta todavía como
cierre definitivo**. Si la parte resuelta también aporta amortiguamiento
apreciable, los datos exigirán ajustar un cierre causal o conservar la dinámica
microscópica correspondiente; no se restará una constante distinta a cada
frecuencia para hacer coincidir las curvas.

La potencia media del puerto, el calor radial positivo del candidato y el
balance del circuito son observables separados. La diferencia entre los dos
primeros **no es** un flujo de calor a los reservorios calculado de manera
independiente. Este ensayo lineal no certifica por sí solo la inyección de
calor B.41 a amplitud finita ni la conservación de energía de un transiente
no lineal. Tampoco fija la cascada fotónica ni las tasas materiales pendientes
de 3.5. Esos límites no se eliminan mediante una etiqueta de etapa cerrada.

## Resultados que deben conservarse

El destino contiene `references/`, `responses/`, los planes exactos y
`workflow_result.json`. Cada respuesta guarda campos complejos, matrices,
momentos energéticos integrados, residuos, índices por energía y el circuito.
Un fallo detiene únicamente su caso y sus trabajadores; los casos independientes
continúan. Una referencia que no converge bloquea solamente sus dependientes.
Las sumas parciales son evidencia de avance, **no un reinicio automático**.

Cuando termine basta avisar en el chat. No hace falta copiar todo el registro.
La interpretación, las figuras comparativas y la decisión de cierre se harán
sobre estos resultados; el texto `COMPLETED` sólo significa que finalizó la
ejecución.

Las ecuaciones y convenciones auditables están en
[cinética](kinetic_derivation.md), [energía y circuito](energy_circuit_derivation.md)
y [ejecución paralela](campaign_runtime.md). Esta respuesta a frecuencia finita
es un diagnóstico del acoplamiento; no sustituye el paso Euler heredado del
solver temporal ni modifica producción o la etiqueta v1.0.0.
