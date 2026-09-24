# Geminga: comandos vigentes

Actualizado el 23 de septiembre de 2026. Cuenta jdiaz, sin administrador.

## Etapa 4: corridas terminadas; no hay otro cálculo largo pendiente

Los dos controles enfocados terminaron en 46,09 min. La malla fina reduce el cambio
de disipación Korzh a 2,55 %, pero el cierre conserva sensibilidad física. La nueva
referencia radial Usadel completó 512 problemas espectrales y 6 candidatos en 4,84 s
con 27 trabajadores y un coordinador. Su fuerza discrepa ampliamente del cierre;
primero corresponde corregir éste, no repetir mallas ni lanzar un fotón pesado.
Véase el [informe vigente](implementation/stage4/followup_20260923/README.md).

## Recursos y paralelismo

Inventario real: 16 núcleos físicos, 32 hilos, aproximadamente 123 GiB de RAM, un nodo
NUMA. El máximo del trabajo es 28 CPU lógicas, dejando dos núcleos físicos completos
libres. No son 16 nodos de cómputo. Se revisan afinidad, cuotas visibles y memoria
al iniciar; el presupuesto usa como máximo 90 % de recursos disponibles. Un grupo
compartido distribuye casos y consultas independientes. Cada proceso usa un hilo
BLAS/OMP; no se multiplican grupos de procesos anidados. La reserva de memoria es
conservadora, no un límite impuesto por el sistema operativo. Se registran
recursos, progreso, tiempo transcurrido y ETA; al inicio la ETA puede estar
dominada por el arranque del grupo y se corrige con el avance.

## Comprobaciones ligeras opcionales

Verificar entrega actual y cadena histórica (segundos, sin cálculos físicos):

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/verify_followup.py
```

Pruebas analíticas y presupuesto de recursos (segundos):

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -m unittest discover \
  -s tests -p 'test_stage4*.py' -v
```

Reproducir la referencia radial ya terminada (opcional; 4,84 s medidos, reservar
hasta 1 min para arranque/carga; hasta 28 CPU lógicas y 30 GiB de reserva conservadora,
menos de 30 MB de salida). No modifica producción ni calcula un transiente.
Cambiar el nombre de salida si ya existe; nunca sobrescribe:

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/radial_usadel_reference.py \
  --plan docs/implementation/stage4/followup_20260923/radial_reference_plan.json \
  --output-root tmp/stage4_radial_reference_reproduction --execute
```

Salidas: identity.json registra fuentes/recursos; progress.jsonl conserva avance;
tasks/ guarda modos espectrales; R8_N128.npz,R8_N256.npz,R12_N256.npz contienen
curvas y estimación de cola separada; summary.json concentra comparaciones.

## Horizonte futuro con fotón

Se simulará hasta un cruce confirmado de Vout más un margen, con un techo finito
si no se observa disparo. Sólo se recorta la ventana temporal: mismas ecuaciones,
material, depósito, geometría, malla, bordes, precisión y circuito completo de la
memoria. No se adelanta la recuperación ni se acortan constantes de tiempo. El
estado final permitirá continuar. Los parámetros del criterio se definirán para
el ensayo; no se identificará la falta de cruce con una latencia ficticia. Esta
política todavía no reemplaza el modo heredado que espera un máximo del pulso.

## Historial íntegro

La [libreta previa](implementation/stage4/followup_20260923/previous_delivery_exact/docs/GEMINGA_COMMANDS.md)
conserva los dos comandos completados y los enlaces a entradas anteriores, sin
alterar su evidencia. No hace falta volver a ejecutar las campañas terminadas.

## Política de cómputo

Cálculos razonablemente previstos de más de cinco minutos se entregan al usuario
con propósito, comando copiable en el chat y esta libreta, salidas y recursos.
No se lanzan desde el agente ni mediante trabajos de fondo o sondeos. Los
chequeos inciertos llevan timeout de 240 s; si vence se conservan como incompletos y
se entrega el comando sin reinicios repetidos o subdivisiones para eludir el límite.
No se requieren nuevas sesiones ni comandos screen.
