# Revisión independiente del cierre de la etapa 1

**Dictamen: se cierra la admisión electrónica del catálogo y sus ensayos reducidos de uso. La entrada fonónica queda admitida únicamente como forma numérica condicionada.** El certificado `closure_admission.json` vincula los archivos exactos mediante SHA-256. No autoriza promoción a producción ni declara completado el acoplamiento dinámico íntegro de D.4, punto 2.

Los criterios se fijaron antes de los resultados. Una aclaración de diseño previa a la primera corrida permitió reconstruir linealmente los niveles para el transporte, conservando el criterio de sesgo respecto del catálogo nativo. El contrato inicial queda archivado en `review/acceptance_criteria_initial.json`; las tolerancias no se ajustaron para acomodar los resultados.

## Resultado electrónico

Pasaron las **53 puertas** de `review/electronic_assessment.json`. Se usaron el catálogo R2 y su código originales, verificados por hash y por una repetición nueva de sus comparaciones físicas/API. Las trayectorias se repitieron con 80, 160 y 320 pasos; el transporte empleó 129, 257 y 513 puntos compartidos de energía.

| Comprobación | Resultado en la resolución final |
|---|---:|
| BGK: deriva energética escalada | 2,55 × 10⁻¹⁴ |
| Espectro móvil: residuo con trabajo de campos | 1,53 × 10⁻¹³ |
| Celda con amplitud propia y calor del condensado | 5,55 × 10⁻¹² |
| Transporte: residuo de energía reconstruida | 2,22 × 10⁻¹⁶ |
| Transporte: sesgo de energía frente a R2 nativo | 0,0212% |
| Transporte: error de remapeo FD | 0,0688% |
| Transporte: deriva del número de cuasipartículas | 0,00233% |

Las poblaciones se mantuvieron entre cero y uno sin recortes. BGK aumentó la entropía en los casos no térmicos y conservó su energía; su número de cuasipartículas puede cambiar, como exige el cierre elegido. La fuente deposita la potencia indicada desde el vacío. En la celda aislada, la amplitud pasó de 0,6 a 0,88582 mientras su disipación alimentó las ocupaciones una sola vez. El tiempo y la movilidad son sintéticos.

La auditoría independiente añadió inversión FD por bisección, una trayectoria sin corriente con niveles/fuerzas BCS analíticos e integrador DOP853, y comparación del intercambio de cada par con una exponencial matricial. El control BCS conservó energía a 1,06 × 10⁻¹⁵. En transporte, las identidades de momentos de la reconstrucción se verificaron a 1,43 × 10⁻¹⁴; una cola privada con población no térmica permanece desconectada. Se corrigieron antes del cierre la extrapolación silenciosa, el peso del nodo que termina el soporte común y una cancelación que podía producir una población negativa de redondeo.

La energía conservada por el transporte es la de su **reconstrucción**, no exactamente la cuadratura nativa R2. Su sesgo y el defecto de conteo son resultados separados, ambos menores que 0,1%. El extremo se detiene en el último nodo real y no extrapola hasta el corte nominal: el intervalo no representado y la cota de energía descartada por debajo del umbral inferior se publican en `cells_results.json`. Estos resultados no prueban una cota uniforme para poblaciones arbitrarias.

## Resultado material

La repetición independiente reconstruyó exactamente el CSV derivado y la correspondencia con sus filas de origen. Las integrales cerradas por segmentos de λ y los tres momentos de α²F coincidieron con el cálculo del autor a menos de 4,38 × 10⁻¹⁶ relativo. El cociente entre interpolantes es finito tras imponer soporte común, y α²F = DOS × cociente se verifica a 2,23 × 10⁻¹⁶.

El recorte explícito resulta útil para ensayar una forma espectral trazable, pero **no es un cambio pequeño global** bajo el umbral fijado de 0,1%: λ cambia 0,411%, el momento que normaliza la energía de la burbuja cambia 0,792% y el segundo momento cambia 0,935%. Los perfiles concentrados en la cola o en el hueco interno tampoco permiten concluir inocuidad. La DOS nula con α²F no nulo empieza ya en el eje 7,36881, no sólo en la cola negativa. Si se supone THz, equivale a 30,475 meV y entra en la ventana de reacciones examinada.

La fuente pública sustenta una receta de tratamiento, pero su código espera un archivo de otro nombre que no está publicado. No se ha acreditado la base por átomo/celda ni la normalización absoluta de la tabla numérica disponible. Los resultados en kelvin o meV mantienen explícita la hipótesis de eje THz. La coincidencia aproximada con un conteo esperado de modos o con λ = 1,2 no se usa para fabricar esa evidencia.

## Siguiente implementación admitida

Completar en una y dos celdas las reacciones electrón–fonón con una entrada sintética admisible, acopladas al condensado móvil y al transporte a energía fija. Deben compartir los eventos de intercambio energético y verificar Pauli, positividad fonónica, equilibrio y refinamientos. Esta entrega no evoluciona fonones y no afirma haber probado su positividad. La movilidad material, las unidades fonónicas, el símbolo principal, el núcleo, los reservorios y el circuito permanecen como tareas explícitas antes de los transientes completos. `NEXT_STAGE.md` desarrolla ese orden.
