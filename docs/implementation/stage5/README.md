# Etapa 5: preparación DC interior pendiente; fotón aún no iniciado

La [etapa 4](../stage4/closure_20260924/README.md) se cierra como desarrollo.
Su respuesta reactiva es utilizable como referencia de implementación; el
calor absoluto y la continuidad local del acoplamiento reducido no se promocionan.
No se exige repetir los controles estáticos ni las trayectorias térmicas aceptadas.
Las observaciones posteriores al cierre cambian el próximo objetivo: preparar
una sección interior físicamente admisible del detector antes de inyectar el fotón.
El [contrato DC](prephoton_dc_20260924/physical_scope.md) y su
[registro de fuentes](prephoton_dc_20260924/source_decisions.json) tienen
prioridad para la siguiente campaña. No hay aún resultados de su corrida larga.

## Próximo ensayo: conservar el piso DC de la cinta

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
los 10 ps y las comparaciones completas siguen pendientes. Los cálculos previstos de más de cinco minutos se entregan para
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
