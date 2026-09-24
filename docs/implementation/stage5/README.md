# Etapa 5: entrada preparada, transiente fotónico aún no iniciado

La [etapa 4](../stage4/closure_20260924/README.md) se cierra como desarrollo.
Su respuesta reactiva es utilizable como referencia de implementación; el
calor absoluto y la continuidad local del acoplamiento reducido no se promocionan.
No se exige repetir los controles estáticos ni las trayectorias térmicas aceptadas.

## Primer trabajo: completar la unión física

1. **Continuidad interior.** Usar la misma malla Delaunay-Voronoi y el mismo
   estado de 13,002 microamperios. Reutilizar los operadores dispersos existentes
   para resolver fase y potencial fuera de la base de 6-10 modos, o enriquecerla
   de manera controlada. Comparar corriente compleja integrada en cortes y
   potencia del puerto, además de la señal externa. Un residuo proyectado pequeño
   no sustituye esta prueba. La discrepancia de corte medida es 4,045 % refinada;
   no se atribuye exclusivamente a la base sin un contraste independiente.
2. **Energía y trabajo.** Cerrar el balance independiente de energía interna,
   trabajo de espectro móvil, fuente DC y reservorios y deposición B.41. No definir
   el calor perdido como la diferencia necesaria para que cierre la identidad.
   El factor del calor radial actual fue auditado: dividirlo por dos sería una
   corrección injustificada. Antes de una trayectoria larga, comprobar esta unión
   en una respuesta débil a la misma polarización y circuito.
3. **Transiente débil integrado.** Conectar los bloques admitidos con el Euler
   KWT heredado; evolucionar poblaciones, fase, potencial y circuito conjuntamente.
   Refinar los observables con señal útil, manteniendo el criterio práctico del 2 %
   donde fue declarado. No exigir precisión relativa sobre una cola extinguida.

La etapa no requiere inventar un nuevo integrador ni un ajuste instantáneo del
desequilibrio electrón-hueco. La respuesta armónica no reemplaza este transiente.
Los cálculos previstos de más de cinco minutos se prepararán después del piloto
de coste, con comando en el chat y GEMINGA_COMMANDS.md, barras/ETA y presupuesto
común máximo del 90 % de CPU y RAM disponible.

## Antes del primer fotón de Korzh

- Acreditar la normalización volumétrica de la DOS fonónica NbN y sus tasas.
- Determinar la transferencia tras la cascada y justificar el ancho de la
  preparación gaussiana. No adoptar los intervalos que quedaron sin seleccionar
  en 3.5 ni presentar parámetros ajustados como medidas.
- Registrar material y movilidad coherentes: los controles actuales usan escalas
  D=0,5 cm²/s y R□=608 ohm, pero tiempos KWT heredados de 0,5/2,47 ps, no la pareja
  Korzh 6/24,7 ps. No trasladar sus tiempos de relajación al experimento.
- Mantener inicialmente 2D; justificar cualquier tramo 1D y su transporte a igual
  energía. La caja de 160 x 80 nm es un control, no una longitud fotónica admitida.
- Acordar observable Vout, umbral de cruce, confirmación y margen posterior.
  Preservar circuito completo y constantes de tiempo. La recuperación de
  nanosegundos no es requisito inicial; una corrida sin cruce queda censurada.

La primera comparación experimental sigue siendo latencia relativa de 775 y
1550 nm y formación del hotbelt en la cinta de 80 nm. El reloj del modelo empieza
en la transferencia; el retardo óptico anterior permanece explícitamente desconocido.
