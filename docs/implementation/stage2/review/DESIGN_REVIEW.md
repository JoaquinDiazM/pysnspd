# Revisión previa de la etapa 2

Los criterios se fijan antes de medir las pruebas nuevas. La aceptación se refiere al algoritmo de una y dos celdas; no acredita tasas absolutas de NbN ni el problema espacial de D.36.

## Identidades que deben sobrevivir a la discretización

Cada evento debe compartir una sola tasa entre electrones y fonones. Con capacidades electrónicas `4 w_i`, una dispersión cambia un electrón de estado y una recombinación cambia dos. Para pares no ordenados, el peso diagonal de recombinación es la mitad y su multiplicidad electrónica es dos. En las unidades declaradas, el factor de los pares fuera de la diagonal es `8 pi Delta0 t_ref / hbar`; procede de cambiar las dos medidas energéticas a coordenadas de estados. No se puede comprobar este factor únicamente conservando energía: una tasa multiplicada por dos también conservaría energía.

Si un fonón de energía Ω se reparte entre dos nodos fijos, los pesos deben sumar uno y reproducir Ω. Los productos geométricos de factores Bose conservan la relación de equilibrio porque su cociente es `exp(sum_k b_k Omega_k / T) = exp(Omega/T)`. El mismo argumento se aplica a factores Pauli para un evento de transporte entre estados repartidos. Además, cada evento aislado tiene producción de entropía proporcional a `(A-B) log(A/B)`, no negativa. Estas identidades justifican las pruebas discretas, pero no demuestran que las tasas a resolución finita coincidan con la integral continua.

## Riesgos numéricos concretos

- Las potencias fraccionarias de una ocupación tienen derivadas singulares en cero. El equilibrio interior no prueba el comportamiento desde el vacío: se requieren caras invariantes, trayectorias de borde y refinamiento temporal sin recortar valores negativos.
- El transporte histórico usa otra cuadratura de energía. La etapa 2 debe trabajar con una contabilidad nativa compartida, o mostrar y converger todo defecto de proyección. Sumar dos balances conservativos de energías distintas no crea un balance conjunto correcto.
- La respuesta fonónica de una lista finita de pares es una medida granular. Para comparar rejillas se proyecta conservativamente a bandas físicas comunes; interpolar dos arreglos de derivadas no mide su convergencia.
- Una lista que sólo construye pares dentro del catálogo omite absorciones que terminan fuera. Su contribución debe estimarse desde una referencia independiente, no declararse nula por construcción.
- La cancelación del equilibrio y un residuo energético pequeño pueden coexistir con una tasa equivocada. Se comparan por separado potencias de dispersión/recombinación, respuesta poblacional, transporte continuo, corte infrarrojo y cuadratura electrónica superconductora.

El presupuesto de consistencia espectral es 0,1%; el de error temporal es 0,01%, una década menor. El balance global se exige a 1e-7 en su escala declarada, y las identidades locales a aproximadamente precisión de redondeo. Es una jerarquía numérica para ensayos sintéticos, no una tolerancia de latencia experimental.

## Dictamen del punto de continuación

La etapa 2 queda **PENDING_COMPUTE / NOT_CLOSED**. El dictamen verificable es
`../stage2_admission.json`: enlaza los archivos exactos por SHA y mantiene
separadas las pruebas estáticas, la dinámica y la admisión material. Los
criterios originales y los artefactos de R2 permanecen intactos.

Los controles independientes identificaron que los 180 estados no resolvían
algunas distribuciones, aunque conservaran la energía. La representación
complementaria usa ahora 630 estados y refina a 1260 y 2520. En la resolución
base, el error escalado máximo de los momentos de energía es 2,74·10⁻⁸ y el de
la respuesta a Γ es 4,73·10⁻⁶; ambos están dentro del presupuesto. La comparación
espectral de alta precisión también verifica las derivadas del mismo potencial.
El transporte tiene un error máximo de respuesta de 0,0039463 %, que disminuye
al refinar. La movilidad heredada, las unidades y la entrada Debye sintética
pasan sus controles.

Para colisiones, la configuración propuesta es 630 estados electrónicos,
1025 fonónicos y corte infrarrojo 0,005. El error estático máximo de interpolación
fonónica es 0,038833 %, y el refinamiento interno de la cuadratura de eventos
da 8,01·10⁻⁷ relativo. Son controles específicos: aún falta comparar toda esta
configuración contra la referencia continua independiente al mismo corte.
Las referencias de corte 0,01 tienen errores débiles menores que 0,1 %, pero
una pequeña meseta de BCS no permite afirmar convergencia monótona. Se conserva
una referencia adicional a η finito para distinguir este efecto; cualquier
descomposición sin cerrar permanece pendiente.

Los ensayos temporales terminados con 513 fonones tienen balances pequeños y
poblaciones válidas, pero esa resolución falló el control estático de
distribución. No se trasladan sus resultados a la configuración de 1025 nodos.
La referencia DOP853 bajo límite de 240 s no produjo una trayectoria completa.
La evidencia registra salida vacía y código 1; no identifica una señal concreta
ni permite inferir una inestabilidad física. No se repite ni se fracciona ese
cálculo: la continuación corresponde al comando manual de Geminga.

Quedan pendientes la referencia temporal, tres resoluciones de tiempo, la
convergencia de poblaciones durante la dinámica, los casos de borde/equilibrio
completos y los controles de fuerzas y soporte sobre las trayectorias finales.
Las 509 regresiones pasan en 46,76 s con las fuentes verificadas; esto confirma
el punto de código conservado, sin sustituir esas pruebas físicas.

No hay promoción a producción ni cierre de la etapa espacial. Debye conserva
su carácter sintético y la forma de NbN sigue condicionada por las unidades,
la base y la normalización absoluta pendientes.
