# Primero completar la etapa 2

La etapa 2 permanece pendiente. La referencia manual corta ya terminó y RK4
pasó contra ella. Ejecutar el lote de `resume_20260921/command_addendum.md`.
Faltan el control temporal independiente, los refinamientos dinámicos sobre la
configuración de 630/1025 estados y los límites de soporte en esas trayectorias.
La comparación estática continua al corte 0,005 ya pasó. Conservar
el sondeo incompleto; no convertirlo en aprobado ni reintentarlo por fragmentos.

El informe final y el avance espacial dependen de esas puertas. El plan que
sigue queda preparado para cuando `stage2_admission.json` permita continuar.

## Plan espacial condicionado al cierre

El dictamen de esta entrega se encuentra en `stage2_admission.json`. Ningún
resultado de una o dos celdas acredita por sí solo un transiente completo de
SNSPD. La siguiente etapa requiere una nueva inscripción de casos y tolerancias.

## Entradas que deben conservarse

- El catálogo R2, sus unidades, el potencial de vacío y su evidencia histórica.
- La malla ocupacional complementaria explícita que admita esta entrega. Los
  180 estados anteriores no se consideran suficientes para toda distribución
  cinética sólo porque sus momentos energéticos sean precisos.
- Un único potencial para energía, fuerzas y trabajo espectral. Los operadores
  de colisión y transporte comparten sus cambios de energía; no hay corrección
  posterior del balance ni recorte de ocupaciones.
- Las ocupaciones electrónicas y fonónicas son variables dinámicas. La
  temperatura equivalente se calcula a partir de la energía; no reemplaza las
  distribuciones por un estado térmico.
- La movilidad KWT y el tiempo de relajación cinética son parámetros distintos.
  El calor de disipación del condensado entra una sola vez en D.18.

## Orden de incorporación

1. Construir un funcional espacial discreto cuya variación produzca todas las
   fuerzas del condensado y la corriente. Comprobar D.36, la rigidez espacial,
   y la cancelación del trabajo entre condensado y campo electromagnético antes
   de integrar. El acoplamiento de esta entrega transporta cuasipartículas;
   todavía no representa gradientes de amplitud o fase entre las dos celdas.
2. Añadir las condiciones de borde, interfaces 2D/1D y reservorios de D.1. Fijar
   explícitamente qué energía y partículas intercambia cada borde. Ensayar un
   estado uniforme, un gradiente débil y una interfaz con espectros distintos.
3. Incorporar Poisson, corriente normal y las ecuaciones circuitales D.28-D.33
   usando un registro conjunto de trabajo eléctrico, calor y energía inductiva.
   Repetir primero equilibrio y respuesta pequeña sin fotón.
4. Ensayar un depósito de energía localizado con entradas Debye sintéticas.
   Separar convergencia espacial, temporal y espectral antes de comparar formas
   de pulso. La admisión material de NbN sigue siendo un requisito independiente
   para atribuir tasas o latencias absolutas al dispositivo.

## Condiciones para avanzar

No promover a producción antes de comprobar convergencia espacial, conservación
con bordes/circuito y estabilidad de la formulación espacial. Un cambio del
número de estados, las cuadraturas o el soporte requiere repetir los controles
cinéticos pertinentes. Los casos fallidos de esta entrega se conservan como
regresiones y ejemplos de resolución insuficiente.

Medir coste de un paso antes de lanzar una trayectoria mayor. Registrar en
`/home/jdiaz/GEMINGA_COMMANDS.md` todo cálculo previsto por encima de cinco
minutos y esperar la ejecución del usuario. El tiempo de los ensayos reducidos
no estima por sí solo el coste del sistema espacial completo.
