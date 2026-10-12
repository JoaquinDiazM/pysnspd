# Etapa 5: preparación DC admitida; integración fotónica pendiente

## Decisión tras el RHS fotónico inicial completo del 8 de octubre

El lote de 224 consultas termina correctamente en 45,8 s, pero no admite cerrar la etapa 5 ni su publicación como modelo concluido. El trabajo K1 es el 92,59/89,40 % del K0 para 775/1550 nm; el defecto independiente es 14,271/8,725 nW. Las derivadas Moyal dominan la corrección, y faltan neutralidad/gauge, los puertos de corriente y el balance temporal conjunto. La [decisión vigente](completed_photon_rhs_20261008/closure_decision.json) conserva las hipótesis de material/fuente y márgenes 1 %/2 %. No repetir DC, UV ni los controles del actuador; no corregir con clipping, corriente óhmica añadida o calor igual al residuo.

La siguiente iteración es de implementación: cierre microscópico inicial G0+G1=0, Gamma_res=0 como candidato explícito, con KWT sólo como precondicionador numérico y no una segunda disipación física. No se presume bien condicionado ni admitido; requiere cortes superconductores dinámicos, neutralidad y puerto CM del mismo sistema, energía al mismo orden y fuentes consistentes. No hay un photon-run del nuevo backend que pueda ejecutarse rigurosamente todavía. Las pipelines nuevas revisan y plotean los datos iniciales guardados, con cero pasos físicos. Se conserva provisionalmente la ventana de 200 ps sin gatillo y 5 ps extra, umbral positivo 0,1 mV, sin reset.

El [informe de esta revisión](../../../output/pdf/implementation/Informe_revision_RHS_fotonico_completo_20261008.pdf) reúne cuatro páginas de plots físicos y la lectura SNSPD. Las pipelines `pipelines/04_initial_photon_rhs_review.py` y `plot_pipelines/04_plot_initial_photon_rhs_review.py` verificaron los datos guardados y producen diez figuras, sin evolución temporal; el [recibo](completed_photon_rhs_20261008/pipeline_review_receipt.json) registra unidades, comprobaciones y recursos. No hay cálculo pesado pendiente del usuario.

## Decisión anterior tras las cinco campañas (histórica)

Se cierran los controles estáticos y de transporte examinados. Ambas repeticiones llegan a 2 ps; las tres nuevas cuadraturas dan un defecto estable de 1,70426/1,64592 pW. No repetir Gauss, UV ni los holds DC admitidos. La cola Matsubara analítica diagnostica el corte térmico, pero completar también la fuerza cambia la acción: no se presenta como reparación del RHS original. Para la siguiente exploración se conserva la fuerza finita F256 de la referencia DC.

La etapa 5 completa **sigue abierta**: falta integrar gap, ambas distribuciones, neutralidad, fuentes físicas, fonones y CM en un mismo RHS. Se prepara un único control inicial real para 775/1550 nm, con fuerza K0 y movilidad KWT heredada, y K1/carga/corriente/trabajo como diagnósticos. No es una trayectoria ni una hotbelt. El límite práctico futuro es 1 % de energía retenida en balance integrado y 2 % de cambio en observables refinados, conservando soporte físico y dirección DC. No se añade calor igual al residuo.

La decisión trazable está en `final_review_20261008/closure_decision.json`; el análisis entregable es `output/pdf/implementation/Informe_final_etapa5_20261008.pdf`, junto a tres documentos matemáticos de ecuaciones principales, auxiliares y precomputadas. Todos los datos pesados nuevos y renders van a scratch; ambas carpetas tmp son enlaces a ese disco.


El piloto inicial de 775 nm midió ||R1||F/||R0||F=1,218 cerca del gap (máximo nodal 4,469). El lote posterior de 224 consultas ya terminó y su revisión vigente está al comienzo de este documento; no repetir ese comando. La razón local no se interpreta como error de latencia y no hubo avance temporal.

## Revisión anterior del 8 de octubre (histórica; sustituida por la decisión superior)

Las cinco campañas preparadas el 5 de octubre ya fueron ejecutadas por el usuario y revisadas. No repetir DC ni el corte UV. La etapa 5 completa permanece abierta porque el balance a energía finita cambia con la cuadratura cerca del borde espectral. El transporte h5 pasó a ambos pasos; las otras dos geometrías requieren completar el ensayo con el nuevo paso común que preserva soporte Pauli sin clipping. Los comandos realmente necesarios están en `/home/jdiaz/GEMINGA_COMMANDS_PRESENT.md`; el historial conserva las entradas usadas.

Se implementaron momentos SI completos de Nambu, su contraste estacionario directo con los espectros guardados, una corrección causal explícita del regulador numérico, archivo de snapshots/checkpoints y preparación de fuente óptica conservativa. La corrección causal no es calor material y sus trazas L/T son cero; no elimina el defecto de trabajo medido. El control estacionario evalúa neutralidad, fuerza y corriente independientes y el circuito CM sobre el estado original; no es una trayectoria nueva acoplada.

El usuario autorizó **material efectivo acústico para pruebas exploratorias** y **ancho gaussiano de referencia con sensibilidad**. La referencia es sigma=7,549834435 nm, derivada de D=5e-5m²/s y una escala de cascada de 0,57 ps de Vodolazov2017, con limitaciones explícitas para 775/1550 nm. Material, depósito y observable quedan en [el contrato físico fechado](photon_readiness_20261008/source_and_observation_contract.md). El espesor es 7 nm; h=5/3,5 nm sólo cambia la malla.

El techo sin gatillo confirmado sigue siendo 200 ps. La observación exploratoria usa 0,1 mV positivo, confirmación 0,5 ps y margen 5 ps; conserva los tres estados CM y aborta antes de aceptar inversión global de corriente. No se admite por estos controles el backend fotónico no térmico general, ni se prepara un reset eléctrico completo. Los mapas DC/SS y la fuente inicial se distinguen de una futura trayectoria real de formación de hotbelt. La entrega de análisis es un nuevo PDF compilado con magnitudes físicas, preservando el informe anterior.


La [etapa 4](../stage4/closure_20260924/README.md) se cierra como desarrollo.
Su respuesta reactiva es utilizable como referencia de implementación; el
calor absoluto y la continuidad local del acoplamiento reducido no se promocionan.
No se exige repetir los controles estáticos ni las trayectorias térmicas aceptadas.
Las observaciones posteriores al cierre cambian el próximo objetivo: preparar
una sección interior físicamente admisible del detector antes de inyectar el fotón.
El [contrato DC](prephoton_dc_20260924/physical_scope.md) y su
[registro de fuentes](prephoton_dc_20260924/source_decisions.json) tienen
prioridad para la siguiente campaña. La corrida larga terminó: la [revisión actual](dc_admission_20261005/README.md)
revalida `DC_PREPARATION_ADMITTED`, cuatro referencias y cinco holds de 10 ps.
No repetir esta campaña. La etapa 5 completa sigue abierta.

## Ensayo completado: conservar el piso DC de la cinta

1. **Equilibrio intrínseco con corriente.** Obtener el gap y espectro homogéneos
   a la corriente elegida y usar ese mismo estado en los cortes longitudinales.
   La magnitud de $\Delta$ debe ser un piso uniforme; su fase varía como $qx$.
   No usar terminales metálicos ni espectros de corriente cero para esta ventana
   interior. Comparar longitudes a igual corriente, con avance de fase $qL$.
2. **Relajación sobre la malla dual.** Preparar la referencia autoconsistente
   con [prepare_dc_reference.py](../../../sandbox/stage5_prephoton/prepare_dc_reference.py),
   conservando los operadores Delaunay-Voronoi. Medir el perfil de amplitud,
   la corriente en cortes interiores y la influencia de longitud y resolución.
   Una solución analítica copiada a los nodos no acredita equilibrio discreto.
3. **Conservación temporal sin fotón.** Usar el paso KWT Euler heredado y el
   circuito CM en [dc_hold.py](../../../sandbox/stage5_prephoton/dc_hold.py),
   con fuente DC y sin una sonda AC adicional. El campo normal de este control
   usa el cierre óhmico heredado y las poblaciones se fijan en equilibrio:
   no se presenta como evolución cinética general fuera del equilibrio.
   Medir la deriva de amplitud, corriente y señal basal, además del error
   de conservación. La fase de puerto puede evolucionar coherentemente con
   su potencial; no se suprime la respuesta del circuito fijando ambos extremos.

El [corredor conjunto](../../../sandbox/stage5_prephoton/run_prephoton_dc.py) y
el [plan de ejecución](prephoton_dc_20260924/README.md) están preparados:
primero prepara las referencias y después ejecuta sólo las conservaciones que
dependan de referencias admitidas. La preparación no es tiempo físico; sólo el
segundo bloque integra el tiempo. El plan y el comando definitivos acompañan
la campaña. El piloto estático y un piloto temporal de 50 pasos ya pasaron;
los 10 ps y las comparaciones completas ya pasaron. La
[continuación vigente](continuation_20261005/README.md) prepara las interfaces
que todavía faltan antes del primer fotón. Los cálculos previstos de más de cinco minutos se entregan para
ejecución manual, con comando en el chat y GEMINGA_COMMANDS.md, barras/ETA y
presupuesto común máximo del 90 % de CPU y RAM disponible.

## Mismo circuito, inductancia exterior fija

Se mantienen las tres variables $(I_b,I_s,v_c)$ y la topología CM. La fuente
constante se inicializa como $V_b=R_bI_{\rm DC}$ y la señal basal es cero.
Para el escenario motivado por Korzh se fija

$$L_{k,\rm ext}^{\rm ref}=96\,\mathrm{nH}
+(5\,\mathrm{\mu m}-L_{\rm res})\mathcal L_{\rm bulk}^{\rm ref},
\qquad \mathcal L_{\rm bulk}^{\rm ref}
=\frac{\hbar}{2e}\left.\frac{dq}{dI}\right|_{I_{\rm DC}}.$$

La derivada se evalúa una sola vez por estado inicial. No hay $L_k(t)$ ni
repartición instantánea de inductancia. Los 96 nH pertenecen al inductor añadido
exterior; no se les resta otra vez el segmento resuelto. La parte activa no
resuelta se contabiliza de forma compatible con el material elegido. Los
100 nH aproximados obtenidos de 64 pH/cuadro son una comparación de escala,
no una medición del total ni una segunda inductancia que se suma al solver.

La ventana representa un interior recto del tramo activo de 5 µm. Se fija la
coordenada longitudinal del futuro fotón y se excluye el jitter geométrico
longitudinal. El estudio posicional posterior corresponde al ancho de 80 nm.
El resto de los parámetros CM conserva su procedencia de la memoria; no se
afirma haber reconstruido toda la electrónica experimental de Korzh.

## Antes del primer fotón de Korzh

- Completar el balance independiente de energía interna, trabajo espectral,
  fuente DC, reservorios y deposición B.41 para una perturbación finita. Un
  equilibrio sin disipación no corrige por sí mismo la discrepancia del calor
  encontrada en etapa 4. No se define el calor perdido como el residuo necesario
  para hacer cerrar la identidad.
- Conectar evolución de poblaciones y transporte a los puertos del circuito
  para el transiente fotónico. El control térmico DC no elimina el desequilibrio
  electrón-hueco del modelo general ni admite su relajación instantánea.
- Acreditar la normalización volumétrica de la DOS fonónica NbN y sus tasas.
- Determinar la transferencia tras la cascada y justificar el ancho de la
  preparación gaussiana. No adoptar los intervalos que quedaron sin seleccionar
  en 3.5 ni presentar parámetros ajustados como medidas.
- Registrar material y movilidad coherentes: los controles de etapa 4 usan escalas
  D=0,5 cm²/s y R□=608 ohm, pero tiempos KWT heredados de 0,5/2,47 ps, no la pareja
  Korzh 6/24,7 ps. No trasladar sus tiempos de relajación al experimento.
- Mantener inicialmente 2D; justificar cualquier tramo 1D y su transporte a igual
  energía. La caja de 160 x 80 nm es un control, no una longitud fotónica admitida.
- Acordar observable Vout, umbral de cruce, confirmación y margen posterior.
  Preservar circuito completo y constantes de tiempo. La recuperación de
  nanosegundos no es requisito inicial; una corrida sin cruce queda censurada.

Las respuestas armónicas de etapa 4 se conservan como evidencia histórica; no
se exige una nueva campaña AC ni otras perturbaciones ajenas al siguiente
problema físico. La primera comparación experimental sigue siendo latencia relativa de 775 y
1550 nm y formación del hotbelt en la cinta de 80 nm. El reloj del modelo empieza
en la transferencia; el retardo óptico anterior permanece explícitamente desconocido.
