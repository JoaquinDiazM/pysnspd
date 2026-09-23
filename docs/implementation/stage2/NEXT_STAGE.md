# Continuación tras el cierre de desarrollo de la etapa 2

**22 de septiembre de 2026. Etapa 3: PREPARADA_NO_INICIADA.**

La etapa 2 se cierra en su **alcance de desarrollo**, con la evidencia disponible
y según la [decisión de cierre](closure_20260922/closure_decision.json).
Su certificado estricto de convergencia dinámica
entre mallas sigue incompleto: el cierre no cambia los resultados originales a
PASS ni equivale a admisión material o de producción.

La [revisión de resultados guardados](practical_review_20260922/saved_results_audit.json)
recoge 21 tareas completadas y 13 trayectorias del método protegido. Sus
comparaciones temporales propias de una y dos celdas, balances y poblaciones
pasaron en los casos registrados. El
[par fonónico de 2049 nodos](practical_review_20260922/raw/two_guarded_ph2049_pair_320_640.json)
conserva `FAIL_PAIR_DIAGNOSTIC`: diferencia `4.565195781e-5` frente al presupuesto
auxiliar `2.5e-5`. No se infiere convergencia de las mallas pendientes y no se
ordena relanzar automáticamente el lote anterior.

El [dictamen de etapa 2](stage2_admission.json) distingue ese alcance del
certificado numérico incompleto. `numerical_admission = false` **ya no funciona
como bloqueo global de la siguiente etapa**. Cada ensayo futuro deberá fijar
sus observables y controles pertinentes antes de correr; las incertidumbres
heredadas limitan sus conclusiones y no se ocultan.

La dinámica excitada anterior mantuvo `Gamma = 0`, y el equilibrio de dos
celdas, `Gamma = (0, 0.1)` fijo. No verificó fase, potencial, energía espacial
ni circuito. La etapa 3 deberá comprobar por sí misma el trabajo espectral
con `Gamma(q_delta)` variable y el trabajo eléctrico de puerto usando fuerzas
y corriente de la misma energía; no se transfiere el balance de amplitud sola.

## Siguiente implementación, todavía no iniciada

El [plan de etapa 3](../stage3/README.md) y su
[contrato preparado](../stage3/entry_contract.json) desarrollan D.4.3:
**resolver el ensayo espacial admitido con bordes y circuito**. El orden es:

1. **Funcional estático:** energía discreta, fuerzas cartesianas y corriente de
   la misma energía; soporte y estabilidad D.36. Primer piloto 1D con poblaciones
   congeladas y cierre periódico de diagnóstico.
2. **Bordes y reservorios:** sustituir ese cierre por paredes aislantes,
   interfaces 2D–1D, continuaciones y reservorios; comprobar flujos y su trabajo.
3. **Carga y potencial:** resolver conservación de corriente y potencial con
   signos de terminales explícitos, sin introducir una acumulación de carga que
   el modelo continuo no contiene.
4. **Circuito de la memoria:** acoplar las tres ODE de `(Ib, Is, vc)` de la
   [adenda vigente](../../modelo_v0_4/actualizaciones/circuito_memoria_20260922.md),
   con capacitor de lectura, convención pasiva y sólo la inductancia exterior
   no representada ya por los campos. El balance completo es CM.9, que sustituye
   al bloque circuital y balance D.28–D.29/D.33 históricos.
5. **Dinámica débil y depósito local:** primero sin fotón, luego con una entrada
   sintética localizada dentro del dominio admitido. Medir la precisión en los
   observables elegidos y distinguir los errores espacial, temporal y espectral.

El ensayo final tendrá terminales, reservorios y circuito. La periodicidad de
3A sólo aísla el funcional; no representa los bordes del detector. Se conserva
el control negativo D.36 (`|Delta|/Delta0=0.60`, `q*ell0=1`, `delta/Delta0=0.10`)
y se rechaza ese estado para evolución, sin recortar su rigidez.

## Entradas y criterios de la continuación

Se mantienen catálogo R2, procedencia, energía común y representación
ocupacional explícita. Las interfaces transportan a energía física común y
contabilizan el trabajo espectral; KWT y relajación cinética siguen siendo
cierres distintos, y `Q_Delta` entra una sola vez. Para corriente positiva de
izquierda a derecha, `Vdev=phi_L-phi_R` e `Is*Vdev` es trabajo que entra en el
dominio resuelto. No se añaden los 10 nH **totales** de la memoria sobre la
energía de superflujo ya resuelta.

Geometría, campos, condiciones, normas de error, referencias y presupuestos
se fijarán antes de cada ensayo. Los criterios nuevos serán proporcionales a
lo que se quiera distinguir; no se imponen ahora tolerancias universales ni se
reclasifica retrospectivamente el fallo fonónico. Los problemas estructurales
invalidan el cálculo afectado; la precisión insuficiente se informa y delimita
la conclusión. Una etapa exploratoria no certifica umbrales, latencias ni NbN.

El benchmark de etapa 2 con 320 pasos en tiempo normalizado `0–2` no prescribe
el paso de la futura PDE. Primero se mide el coste y se estima la resolución
necesaria. Todo cálculo previsto de más de cinco minutos queda registrado en
`/home/jdiaz/GEMINGA_COMMANDS.md` y espera ejecución del usuario. No se implementa
ni ejecuta etapa 3 en esta entrega; núcleo, transitorios completos y producción
conservan los alcances posteriores D.4.4–D.4.5.
