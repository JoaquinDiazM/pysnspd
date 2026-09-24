# Ejecución paralela de los próximos casos acoplados

`sandbox/stage4_core/coupled_campaign.py` sólo coordina programas ya definidos.
No contiene ecuaciones del detector, tolerancias físicas ni un solver de
reemplazo. Los argumentos del solver acoplado se incorporarán al plan cuando
su interfaz esté fijada. Este protocolo no solicita todavía una simulación.

## Un presupuesto común comprobado al ejecutar

El runner reutiliza `parallel_runtime.py` para leer topología, afinidad,
cuotas de CPU y memoria disponible. Asigna como máximo el 90 % de esos
recursos y reserva núcleos físicos completos. En la topología de Geminga de
16 núcleos/32 hilos, el techo es 28 CPU lógicas, incluidos todos los
coordinadores. BLAS, OpenMP y demás bibliotecas numéricas quedan en un hilo.

Cada caso recibe CPU distintas. Su asignación incluye un coordinador y el
número de trabajadores que sustituye `{workers}` en el argv. Dos casos que
ocupen todo el presupuesto pueden usar 12 y 13 trabajadores, sus dos
coordinadores y el coordinador de campaña: 28 en total. Un caso aislado puede
usar 26 trabajadores y dos coordinadores: también 28. Es un techo; el plan
puede limitar los trabajadores y dejar más capacidad libre.

La reserva conservadora de RAM es 3 GiB para el coordinador de campaña,
3 GiB para cada coordinador de caso y 1 GiB por trabajador. Si no cabe en el
90 % disponible, disminuyen los trabajadores o los casos simultáneos. Es una
reserva de planificación, no un límite de memoria impuesto por el sistema.
El solver hijo debe respetar `{workers}` y la afinidad heredada; no puede
iniciar por su cuenta otro grupo de trabajadores fuera de esta cuenta.

## Contrato del plan

El JSON usa `schema: "pysnspd.stage4.coupled_campaign.v1"` y contiene:

| Campo | Significado |
|---|---|
| `maximum_parallel_cases` | Máximo de casos independientes simultáneos. |
| `maximum_workers_per_case` | Techo opcional de trabajadores numéricos por caso. |
| `heartbeat_seconds` | Intervalo de impresión, entre 0,1 y 60 s; por defecto 5 s. |
| `working_directory` | Directorio existente de los programas; por defecto el repositorio. |
| `cases` | Lista de objetos con `id`, `argv`, `estimated_seconds` opcional y `expected_outputs` opcional. |

`argv` es una lista exacta de argumentos, nunca texto de shell. Se sustituyen
únicamente `{python}`, `{repo}`, `{workers}` y `{output}`. Los dos últimos son
obligatorios para hacer explícitos el presupuesto y la salida propia de cada
caso. Los caracteres de shell dentro de un argumento permanecen literales.
No se ejecuta `bash -c`, ni se agregan comillas o escapes a los argumentos.

El directorio `{output}` será `<output-root>/cases/<id>` y **no se crea antes
de llamar al solver**. Así puede usarse con programas que rechazan salidas
existentes. Los nombres de caso son únicos y no contienen rutas. Cada entrada
de `expected_outputs` es un archivo relativo a ese directorio; un retorno
cero sin esos archivos queda registrado como `MISSING_OUTPUT`.

La inspección predeterminada no lanza casos ni crea directorios:

```bash
python -m sandbox.stage4_core.coupled_campaign --plan PLAN.json --output-root OUTPUT --dry-run
```

La ejecución requiere añadir `--execute`. El comando concreto se entregará
cuando exista el plan físico; si se estima que supera cinco minutos, debe
aparecer tanto en el chat como en `/home/jdiaz/GEMINGA_COMMANDS.md`, y lo
ejecutará el usuario. El runner no cambia esa política.

## Progreso, ETA y evidencia

La terminal imprime una barra y eventos JSON. El ETA utiliza los tiempos
restantes reportados por los casos o sus `estimated_seconds`; es una estimación,
no una garantía. Antes de disponer de una base se muestra `unknown`. Los casos
pueden imprimir objetos JSON con `fraction` o `progress_fraction` en [0,1],
`eta_seconds`, o la pareja `time_ps` y `horizon_ps`. El resto de su salida se
preserva igualmente en el archivo de consola.

La campaña guarda `executed_plan.json`, `campaign.json`, `progress.jsonl`,
`case_records/<id>.json` y `logs/<id>.console.log`. Cada registro conserva argv,
afinidad, trabajadores, tiempos, retorno, estado y hashes de salidas declaradas.
Los manifiestos se escriben mediante reemplazo atómico después de vaciar a
disco el archivo temporal. Las salidas del solver permanecen bajo `cases/`.

Un caso fallido se intenta una sola vez y no impide que terminen los demás.
El resultado global será `COMPLETED_WITH_FAILURES`; el retorno de campaña es
1. Una campaña íntegramente correcta retorna 0. `SUCCEEDED` significa retorno
cero y presencia de los archivos declarados: la aceptación física permanece
en el resultado de cada solver y no la deduce este coordinador. Un directorio de campaña
existente se rechaza: no hay sobreescritura, reanudación o reintento silencioso.
Al interrumpir la campaña se terminan los grupos de procesos activos, incluidos
sus trabajadores, y se conserva el manifiesto como `INTERRUPTED`.

Los tests usan programas Python diminutos que duermen fracciones de segundo;
comprueban afinidad real, paralelismo, fracaso de un hermano, retorno, ausencia
de resultados, preservación de argv y rechazo de reutilización. No resuelven
física ni suponen un aumento de velocidad paralelo no medido.
# Admission inherited by a nested case

The coordinator passes `PYSNSPD_SHARED_ALLOCATION`, a JSON envelope with its
PID, the case CPU set and worker allowance, and the original shared budget.
The response and stationary-reference runners verify that envelope against
their actual affinity, parent PID, current quota and available RAM. They then
pin their coordinator and each numerical worker separately. The original 90%
reservation is applied once to the whole campaign; a nested case does not
reserve another 10% of an already admitted subset. A standalone invocation
performs its own fresh admission.

`coupled_response.py` writes the final status into `manifest.json`, records
source and result hashes, and saves each worker's CPU affinity, CPU time and
peak resident memory. A failed query is retained alongside successful sibling
queries. `COMPLETED_WITH_FAILURES` returns exit code 1. These are execution
statuses, not physical acceptance decisions.

`biased_strip_reference.py` uses the same runtime protocol for stationary
Matsubara queries. Its sweep count is an optimization index, never physical
time. A consistent last `reference.npz` is always distinguished from the
separate `next_delta_bar` proposal. Reaching the iteration budget without the
declared gap criterion produces `INCOMPLETE_MAXIMUM_ITERATIONS` and exit code
2; it cannot serve as an admitted reference merely because a file exists.
