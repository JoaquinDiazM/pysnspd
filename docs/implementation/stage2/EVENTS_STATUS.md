# Eventos electrón–fonón: evidencia y pendientes

**La etapa 2 no está cerrada.** Los controles estáticos siguientes no sustituyen
la referencia temporal que quedó incompleta ni autorizan producción o una
calibración absoluta de NbN. Los criterios originales no se modificaron.

El candidato conserva energía y número con los mismos pesos que usa en sus
actividades de Pauli y Bose. Su cuadratura resuelve tanto los nodos electrónicos
como los intervalos de energía fonónica. Para Γ=0, los estados virtuales usan
las expresiones causales BCS exactas; para Γ>0 se conserva una reconstrucción
lineal cuyo error se midió en el caso sin gap. No se recortan poblaciones ni
coeficientes, y los eventos acoplados sin soporte se rechazan.

| Control | Resultado medido | Alcance |
|---|---:|---|
| RHS electrónico frente a referencia independiente, 630 estados | 2,994×10⁻⁴ | Tres campos y tres perfiles; Ωmín=0,01 |
| Mismo control, 1260 / 2520 estados | 2,828×10⁻⁴ / 2,915×10⁻⁴ | Hay un piso BCS frente a η=0; no se declara convergencia de ese piso |
| Fonones Lobatto, 513 nodos, 33 funciones de prueba | 1,278×10⁻³: **falla** | Ωmín=0,005; conservado como antecedente |
| Fonones Lobatto, 1025 nodos, mismo control | 3,883×10⁻⁴ | Cumple el umbral fijo de 10⁻³ |
| Error del RHS electrónico por interpolación Bose, 1025 nodos | 5,031×10⁻⁵ | Misma cuadratura electrónica en ambos lados |
| Cuadratura interna del candidato, orden 2 frente a 4 | 8,007×10⁻⁷ | 630 estados, Lobatto1025, Ωmín=0,005 |
| Pruebas del módulo de eventos | 97 aprobadas | Incluyen soporte móvil, nodos separados por una ULP y caras físicas |

El candidato de 630 estados requiere aproximadamente 0,08–0,10 s para construir
una red y 0,042–0,045 s por evaluación del RHS estático en Geminga. Estos valores
no predicen por sí solos un transiente: la red cambia con el condensado y la
referencia temporal adaptativa puede requerir muchas evaluaciones. La matriz
estática completa de tres mallas tardó 229,70 s y alcanzó unos 4,5 GB de memoria.
No debe repetirse automáticamente con una configuración mayor. La red de
2520 estados tiene alrededor de 14 millones de eventos.

Los antecedentes fallidos se conservan. La cuadratura original por pares de
180 estados falló el control continuo; la reconstrucción lineal de estados
virtuales BCS produjo un error del RHS de 2,766×10⁻³; y los intervalos que no
coincidían con las funciones de prueba electrónicas dieron errores de
autoconvergencia del 14,1 %. Son problemas diferentes del error en el momento
conjugado a Γ que se corrigió refinando la malla de conteo cerca del borde.

La matriz electrónica independiente usa Ωmín=0,01. Para el candidato dinámico
Ωmín=0,005 se midieron la interpolación fonónica y la autoconvergencia interna,
pero su comparación independiente al mismo corte queda explícitamente
pendiente. `selected_event_measurements.json` guarda todos los valores para
hacerla sin repetir la medición. `event_summary.json` reúne máximos, costes y
hashes; los JSON completos conservan las filas por campo, perfil y canal.

Comandos reproducibles ligeros, desde `/home/jdiaz/pysnspd` y con un hilo:

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -m pytest tests/test_experimental_kinetic_events.py -q

env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/check_projected_phonons.py \
  --nodes 129 257 513 1025 --infrared .005 --output tmp/phonon_check.json

/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/summarize_events.py
```

El cálculo temporal pendiente y su política de ejecución se registran en
`/home/jdiaz/GEMINGA_COMMANDS.md`; no deben sustituirse por otro ensayo corto ni
por una tolerancia relajada.
