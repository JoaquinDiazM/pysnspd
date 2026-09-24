# Proyección estática del modo de carga sobre un potencial

Se ensayó `hT(E,i)=chi(E)*v_i` con los 72 espectros y las respuestas cinéticas
ya calculadas. Aquí `chi=∂E tanh(E/(2Tb))`, `Tb=0,9/8,65` y la energía se mide
en `kBTc`. El campo real `v` representa un desplazamiento electroquímico
adimensional: no se identifica todavía su signo con el potencial de la memoria.
No hubo nuevas soluciones espectrales ni evolución temporal.

En cada uno de los seis pares caso/eta, la ecuación resuelta fue

`[Σ w chi L_TT] v = -Σ w L_TL hL`,

con `v=0` en los mismos contactos. Se conservaron el operador cinético, las
sondas radiales y angulares y la cuadratura trapezoidal del ensayo precedente.
La respuesta completa `hT(E,i)` sirve como referencia sobre ese mismo espectro.
Se integraron la fuerza generalizada, la corriente de carga `2Σw*jT` y el flujo
de energía `Σw*E*jL`; este último contiene expresamente el peso energético.

## Resultado

Los seis diagnósticos terminaron en **2,198 s**, con seis trabajadores y un
coordinador, dentro del presupuesto compartido. Los residuos de carga integrados
en nodos libres quedaron por debajo de **1,38×10⁻¹⁶** y la identidad conjunta de
fuerza y corriente por debajo de **2,90×10⁻¹⁷**. Esto acredita la proyección y la
contabilidad algebraica, no la adecuación física de un cierre de potencial único.

| Sonda angular: proyección frente a respuesta completa | Diferencia relativa observada |
|---|---:|
| Corriente de carga integrada | 19,92–20,57 % |
| Fuerza generalizada, densidad nodal | 3,109–3,172 % |
| Flujo de energía con peso E | 0,415–0,430 % |

Las normas de enlace son `sqrt(Σ J²/c)` en caras con punto medio `r≤4ℓ0`;
las normas nodales utilizan las áreas duales. Los mapas guardan también las
componentes de amplitud y torque de fase de la fuerza. Los cocientes no comparan
máximos nodales con normas espaciales.

En la sonda angular, admitir la respuesta de carga completa cambia la corriente
respecto de `hT=0` en una norma equivalente al 50,02–50,34 % de la corriente
completa. La discrepancia de la proyección equivale al 39,82–40,86 % de ese
incremento. En la sonda radial simétrica la respuesta de carga es muy pequeña;
la asimetría la vuelve visible. No basta usar solamente el caso radial simétrico
para decidir la simplificación.

## Límite que obliga a resolver mejor la energía

La misma cuadratura entrega `Σw*chi=1,180225216`, mientras que la integral
analítica en la ventana es `tanh(Emax/(2Tb))=1`: un exceso de **18,02 %**.
No se renormalizó `chi`. El defecto expone que estos 12 puntos no resuelven bien
la forma estrecha de la susceptibilidad térmica. No es lícito descontar ese
porcentaje de la discrepancia de corriente ni atribuir toda la discrepancia
al cierre físico: ambas afirmaciones exigirían datos adicionales.

La siguiente prueba debe comparar mallas energéticas anidadas y desplazamientos
del contorno con los mismos observables. La conservación integrada por sí sola
no acredita la ley óhmica de la memoria; tampoco estos datos gruesos justifican
descartarla. El diagnóstico es estático, lineal y con espectro congelado. No
admite una población finita, relajación instantánea del modo de carga ni una
dinámica electrostática nueva.

## Archivos

- `identity.json`: fuentes exactas, hashes de entradas y presupuesto de recursos.
- `summary.json`: normas, cocientes, residuos y diagnósticos por energía.
- `fields/`: seis mapas con `v` y los momentos integrados de las respuestas con
  `hT=0`, completa y proyectada, todos sobre idéntica cuadratura.
- `charge_effect_comparison.json`: comparación de la discrepancia con el
  incremento producido por admitir el modo de carga completo.
- `verification_receipt.json`: comprobación local de fuentes y seis mapas.

Reproducción con el comando `next_potential_projection.py` y rutas de los datos
congelados. Su función `worker(job)` admite una lista ordenada de energías para
el adaptador de la próxima campaña; el `main` preserva el contrato de 12 puntos
de esta ejecución. Los resultados nuevos requieren una carpeta nueva.
