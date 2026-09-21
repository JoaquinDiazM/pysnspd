# NbN: recorte trazable y límites de uso

El recorte autorizado permite conservar una **entrada derivada para investigación**, con las operaciones y los hashes registrados. El signo negativo del archivo bruto deja de ser un impedimento absoluto. La revisión nueva sí encuentra dos restricciones distintas: el soporte de α²F no coincide con la DOS, y la normalización absoluta de esta tabla continúa sin acreditarse. No se activa el material en producción.

La tabla y los resultados están en [material_nbn_shape_v1.csv](material_nbn_shape_v1.csv), su [manifiesto](material_nbn_shape_v1.manifest.json) y [material_results.json](material_results.json). Las temperaturas y energías de este diagnóstico usan **condicionalmente** el eje en THz esperado por el código fuente. U y C se integran con las ordenadas originales; no son predicciones absolutas en J/m³ o J/(m³·K).

## Qué acredita la fuente

- La consulta nueva devuelve el mismo repositorio público y los mismos bytes de NbN que R2. El [solver, líneas 92–112](https://github.com/qnngroup/proj-KE-solver/blob/5b6bd747f80016da5ccd51db73c110a8ecc6abf6/solver.m#L92-L112) espera columna 2 = α²F y columna 3 = DOS, elimina duplicados conservando el primero y anula F y α² cuando F es demasiado pequeña. Por tanto, el preprocesado tiene un precedente concreto. Ese código lee `nbn-a2f-ph_2.dat`, ausente en el repositorio público y en el inventario acotado de Geminga; no prueba la convención del archivo disponible.
- [f.m, líneas 3–5](https://github.com/qnngroup/proj-KE-solver/blob/5b6bd747f80016da5ccd51db73c110a8ecc6abf6/f.m#L3-L5) reconstruye el acoplamiento como producto de interpolantes de α² y F. El derivado usa interpolantes lineales de α²F y F y después su cociente, manteniendo exactamente la identidad entre ambos. La comparación numérica incluye ambas elecciones, sin afirmar que sean idénticas entre nodos.
- La [figura 2b y el apéndice VII.4 del artículo](https://arxiv.org/html/2501.13791v3#S7.SS4) muestran el sector acústico en meV y DOS en estados/meV, y describen cálculos DFPT. No proporcionan el historial de posprocesado de las tres columnas publicadas. La fuente no acredita que las ordenadas negativas provengan de digitalización, interpolación o extrapolación. El ajuste usado para evitar frecuencias imaginarias no explica por sí mismo ordenadas negativas de la DOS a frecuencia positiva.
- [Quantum ESPRESSO](https://www.quantum-espresso.org/Doc/INPUT_MATDYN.html) documenta DOS en estados/cm⁻¹ normalizada a 3·nat; el [tutorial de EPW](https://docs.epw-code.org/tutorials/tutorial_04/index.html) muestra su conversión a estados/meV. Esas convenciones admiten más de un posprocesado. Interpretar aquí la tercera columna como estados/meV daría 2.92979 modos brutos o 2.97747 tras recorte, pero la proximidad a 3 no es una prueba y **no se adopta ese factor**. También falta casar la base átomo/celda con la densidad usada por el modelo.
- La integral positiva discreta da λ = 1.215665, cercana al parámetro 1.2 del solver. λ no contiene F y es invariante al cambio lineal de escala del eje: esta comprobación no puede identificar la unidad de la DOS.

## Operaciones y resultados

El derivado v1 ordena establemente, conserva el primero de cuatro pares duplicados, hace F = 0 donde F ≤ 0, hace α²F = 0 en ese mismo soporte y fija α²F(0) = 0 explícitamente para que λ sea integrable. No renormaliza F. El cambio del origen elimina apenas 2.754×10⁻⁸⁹ unidades de área; su importancia es definir el límite matemático. Las elecciones primero/último/promedio quedan como controles, no como tres entradas candidatas.

Se eliminan 361 ordenadas negativas. El área DOS pasa de 0.7084200575 a 0.7199501383. Para las temperaturas ensayadas de 0.5 a 80 K, el mayor cambio de U/C es 0.01033 %; a 160 K C cambia 0.3367 %, y a 1000 K 1.554 %. El criterio previo de efecto pequeño era 0.1 %; no se declara un límite continuo de temperatura a partir de esos puntos.

El soporte común retira 0.6085 % del área α²F: λ cambia −0.4112 %, el momento energético M1 cambia −0.7923 % y M2 cambia −0.9348 %. Estos cambios superan el criterio previo. En B.39 debe recalcularse el normalizador de la burbuja: reutilizar el antiguo conservaría 99.2077 % de la energía; usar consistentemente el nuevo conserva 100 %. Esto normaliza la inyección, no el conteo de modos de F.

El cociente bruto α²F/F no está definido en el hueco que empieza en 7.36881 THz (30.47495 meV) mientras α²F = 0.14577, ni en el borde de 16.3917 THz (67.79062 meV). La interpolación derivada lo vuelve finito; el máximo nodal sigue siendo ≈ 1.05×10⁵ en unidades originales. Finitud numérica no acredita una tasa física correcta.

La banda adicional centrada en 7.45 THz, de ancho 0.04 THz, se añadió **después de descubrir el hueco**, sin reemplazar los perfiles fijados de antemano. Su U no cambia al recortar F, pero su momento ∫E·α²F·n cae 96.03 %. Es una demostración directa de por qué sólo mirar capacidad térmica no basta. La banda 18.5 THz tiene energía bruta negativa; se informa el cambio absoluto y no se declara un error relativo físico pequeño.

La cola negativa queda por encima de la ventana electrónica cerrada aproximada de 32 meV, pero el hueco de 30.47 meV sí la intersecta. Dentro de 0–32 meV los momentos M0/M1/M2 cambian −0.5286/−0.7636/−1.0625 %; dentro de 0–30 meV el corte de soporte no cambia esos momentos. Esto no autoriza ignorar fonones incidentes altos capaces de crear excitaciones fuera del catálogo.

## Decisión

Se conserva v1 como **forma experimental condicionada**, adecuada para estudiar sensibilidad al preprocesado y para los regímenes expresamente comprobados. No se certifica un kernel global equivalente al bruto: los controles no térmicos y los momentos globales muestran cambios relevantes. Para U/C, sólo los puntos ensayados hasta 80 K pasan el criterio 0.1 %; para el diagnóstico normal de intercambio con Tph = 4 K pasan Te = 10 y 40 K, no 160 ni 600 K. No son validaciones de un transitorio completo.

La admisión material absoluta necesita una convención de DOS y base de normalización específica del archivo, una densidad consistente y la convención de tasas. El módulo anterior de admisión y los artefactos R2 permanecen intactos.

## Lectura de figuras

`figures/material_support_review`: DOS fonónica publicada y derivada, cola negativa, corte del acoplamiento en el borde acústico y cociente α²F/F. No es DOS electrónica de Usadel. Los segmentos coincidentes se distinguen mediante línea sólida y discontinua; los huecos no se conectan artificialmente.

`figures/material_impact_review`: (a) cambios térmicos; (b) perfiles no térmicos, con barra verde |ΔU|/∫|F|En y naranja |Δ∫Eα²Fn|/∫Eα²Fn; (c) cambios de λ y momentos globales. En las bandas fijas el ancho es 0.4 THz; el asterisco identifica el diagnóstico adicional de ancho 0.04 THz. La banda alta tiene una referencia U no física, por eso se usa la integral del valor absoluto sólo como escala diagnóstica.
