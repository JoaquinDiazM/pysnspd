# Etapa 4: referencia espacial regular y núcleo autoconsistente

La referencia térmica espacial ya recupera la fuerza de Usadel que el cierre
local anterior no reproducía. El siguiente cálculo busca el núcleo
autoconsistente mediante la **misma energía**, sin asignar tiempo físico a las
iteraciones. La etapa 4 permanece abierta; no se ha cambiado producción.

El [informe con plots](https://github.com/JoaquinDiazM/pysnspd/blob/c80c0f8612d8fd5a2b390938ca5ed4b9d58d3caa/output/pdf/implementation/Informe_etapa_4_energia_espacial.pdf)
presenta resultados y siguiente ensayo. El [análisis reproducible](analysis.md)
incluye limitaciones, corte espectral, mapas e identidades.

## Qué cambió y qué se comprobó

El módulo experimental `pysnspd/experimental/thermal_spatial_usadel.py` resuelve
campos espectrales complejos por frecuencia positiva de Matsubara. Una
parametrización regular mantiene su normalización incluso donde Δ=0. Fuerza
y corriente derivan de la misma energía discreta. No se añade K0, no se ajusta
δ y no se fuerza a la fase espectral a coincidir con la del condensado.

La [deducción y revisión de la memoria](memory_and_sources.md), el
[contrato variacional](variational_contract.md) y las
[opciones de implementación](implementation_options.md) explican por qué es
una representación del funcional térmico espacial, y qué error introduce la
discretización. El esquema conserva volúmenes duales y enlaces de calibre;
el adaptador a la malla Delaunay de producción todavía requiere revisar sus
pesos geométricos, sin importar silenciosamente sus correcciones de respaldo.

| Resultado, mismo perfil y suma finita | Valor |
|---|---:|
| Diferencia de fuerza radial, malla 33² | 2,90 % |
| Diferencia de fuerza radial, malla 65² | 0,720 % |
| Diferencia de fuerza radial, malla 129² | 0,180 % |
| Cambio de fuerza al pasar N=128 a 256, malla 129² | 1,87 % |
| Error relativo de variación energética, amplitud | 9,34 × 10⁻⁸ |
| Error relativo de variación energética, fase | 2,26 × 10⁻⁷ |
| Diferencia máxima de fase espectral en el perfil asimétrico | 18,0° |
| Tiempo de 2.048 problemas espectrales | 60,35 s |

La concordancia de 0,180 % corresponde al contraste espacial **al mismo corte**;
no equivale a un error físico total de ese tamaño. Los dos cortes son sumas
truncadas, sin correcciones añadidas a una derivada aislada. El perfil prescrito
no es aún un condensado autoconsistente ni una solución del detector.

En el centro asimétrico Δ=0, pero el propagador anómalo del modo inferior es
finito, |f₀|=0,0705, y la fuerza cartesiana también. Las correlaciones de pares
pueden persistir por la influencia del entorno. El postproceso incluye este
nodo en las normas complejas y conserva las métricas originales que lo
excluían de la proyección radial.

## Campaña siguiente

El [plan autoconsistente](self_consistent_plan.json) usa cuatro casos:
65²/N128, 65²/N256, 129²/N256 y un núcleo inicialmente asimétrico 65²/N256.
Cada corte se relaja independientemente. Todos comparten el mismo grupo de
procesos y el máximo de 90 % de recursos. La referencia actual admite 27
trabajadores y un coordinador, dejando dos núcleos físicos completos libres.

Se alternan el ajuste del espectro y el mínimo exacto del bloque del
condensado. El descenso de este segundo bloque se demuestra y comprueba
algebraicamente. El objetivo operativo inicial es un residuo RMS ponderado
de 0,1 % dentro del núcleo; se guardan además máximos y observables fuera
de esa región. No se declara automáticamente estacionariedad global.

Las barras informan progreso por barrido y tiempo transcurrido. La ETA estima
el coste del presupuesto de trabajo restante, no el momento garantizado de
convergencia. Hay métricas en cada barrido y checkpoints completos en el
primero, cada cinco y al cerrar un caso o presupuesto. La continuación es
explícita, verifica fuentes/plan/datos y exige una carpeta nueva. Puede repetir
hasta cuatro barridos desde el último checkpoint completo.

El piloto de dos barridos completó 1.792 problemas espectrales en 79,60 s.
La energía bajó en los cuatro casos, pero sus residuos de núcleo quedaron
entre 1,85 % y 2,44 %; el piloto no admite todavía ninguno. La continuación
verificada retomará estos datos y guardará sus resultados en el scratch propio
de jdiaz, con hasta 8,51 GiB proyectados de checkpoints sin compresión.

Limitación operativa: si se interrumpe la escritura de un índice de reanudación,
puede quedar un JSON truncado. El programa se detendrá; habrá que recuperar
manualmente el último índice completo, sin seleccionar otro de forma silenciosa.

El comando vigente y la estimación medida están en
[GEMINGA_COMMANDS.md](../../../GEMINGA_COMMANDS.md) y se entregan en el chat.
No se inicia desde el agente una campaña prevista de más de cinco minutos.

## Lo que todavía no se promueve

Estos ensayos son térmicos, estáticos y con trazas espectrales de borde
prescritas. No establecen una barrera, estabilidad dinámica, KWT, fotón,
latencia ni transporte cinético completo. La representación Matsubara no se
puede insertar en B sustituyendo las distribuciones por una temperatura:
su promoción requiere justificar espectro retardado, conteos, energía y
transporte sin doble contabilización. No se descartan los contratos válidos
de etapas anteriores.

Se mantienen el solver temporal y el circuito de la memoria. La futura
ventana hasta el gatillo Vout más margen sigue cambiando solamente el
horizonte de integración, con el resto del sistema intacto.

## Datos y reproducción

`campaign_plan.json` y `raw/stage4_spatial_energy_reference_20260924/` conservan
identidad, resumen, progreso y doce mapas finales. Los 2.048 espectros completos
(235 MB) permanecen en `/home/jdiaz/pysnspd/tmp/stage4_spatial_energy_reference_20260924`;
`remote_mode_integrity.json` acredita sus hashes sin duplicarlos en GitHub.
El piloto inicial de doce modos se conserva en esta entrega. Los informes
anteriores y sus conclusiones permanecen archivados.
