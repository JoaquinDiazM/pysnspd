"""Archive the complete command notebook and prepare a separate concise draft.

This is an offline documentation transform. It never writes the source notebook,
the active mirror, or a remote file, and it never launches a calculation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
EXPECTED_SOURCE_SHA256 = "67b3556fe57507575568139357204d62d9bd4d74eddd68a26de1d74bcdbb03b9"
ARCHIVE = ROOT / "docs/implementation/stage2/review_20260922/GEMINGA_COMMANDS_before_cleanup.md"
DEFAULT_DRAFT = ROOT / "tmp/stage2_review_20260922/GEMINGA_COMMANDS_clean_draft.md"


DRAFT = """# Geminga: comandos vigentes

Actualizado: 22 de septiembre de 2026. Repositorio: `/home/jdiaz/pysnspd`;
revisión de partida: `main`, commit `6117941`. Cuenta: `jdiaz`, sin administrador.

## Estado y siguiente cálculo

La etapa 2 sigue **sin admisión numérica ni promoción a producción**.
La única continuación vigente es el control temporal de una celda descrito
abajo. No ejecuta dos celdas, refinamientos de malla ni un transiente espacial.

- La referencia manual DOP853 terminó; su copia verificada está en
  `docs/implementation/stage2/resume_20260921/manual_reference/`. No repetirla.
- El lote RK4 `tmp/stage2_validation_20260921_203419/` se detuvo en
  `two_electron_4`; se conservan los casos anteriores completados.
- El lote SSP `tmp/stage2_recovery_20260921_ssp/` se detuvo en `one_ssp_40`
  porque su error de energía excedió el umbral, aproximadamente a los 26 s.
  El limitador no se activó; el siguiente lote comprueba el refinamiento temporal.
  No relanzar el lote anterior.
- El lote RK4 `tmp/stage2_validation_20260922_015109/` confirmó la precisión
  temporal en la malla candidata y repitió el fallo de `two_electron_4`. No repetirlo.

Los bloques de lanzamiento anteriores de etapa 2 quedan **superados: no repetirlos**.
El archivo completo, incluidos comandos de PRE, estados estacionarios, fotones
y diagnósticos anteriores, se conserva en el
[histórico anterior a la limpieza](/home/jdiaz/pysnspd/docs/implementation/stage2/review_20260922/GEMINGA_COMMANDS_before_cleanup.md).
Ese histórico es evidencia; no representa una cola de trabajos pendientes.

## Preparación de la terminal

No es necesario activar Conda cuando se utiliza la ruta explícita del intérprete:

```bash
cd /home/jdiaz/pysnspd
NOTEBOOK_PY=/home/jdiaz/.conda/envs/snspd/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Lectura de resultados existentes

Estos comandos sólo leen archivos. Un registro de error o un lote terminado no
equivale por sí solo a la aceptación del modelo.

Últimos eventos de los lotes revisados:

```bash
tail -n 25 tmp/stage2_validation_20260921_203419/batch_events.jsonl
tail -n 25 tmp/stage2_recovery_20260921_ssp/batch_events.jsonl
tail -n 25 tmp/stage2_validation_20260922_015109/batch_events.jsonl
```

Manifiesto del último lote y resumen de la auditoría actual:

```bash
"$NOTEBOOK_PY" -m json.tool tmp/stage2_validation_20260922_015109/batch_manifest.json
"$NOTEBOOK_PY" sandbox/stage2_cells/review_20260922/show_review_status.py
```

Para comunicar el resultado, basta indicar la carpeta del lote y escribir
«reinicia». Conservar JSON, NPZ, recibos y registros completos, también ante un
fallo. No borrar la carpeta ni volver a ejecutar el lote sobre ella.

## Pruebas ligeras reproducibles

Regresiones pequeñas del integrador y sus evaluadores. No ejecutan un lote de
dinámica ni sustituyen sus criterios de convergencia. Coste esperado: segundos,
un hilo; límite preventivo de cuatro minutos.

```bash
timeout --signal=TERM --kill-after=5s 240s "$NOTEBOOK_PY" -m pytest -q sandbox/stage2_cells/recovery_20260921/test_limited_ssp.py sandbox/stage2_cells/recovery_20260921/test_limited_assessors.py
```

Comprobar la revisión y los archivos modificados antes de comparar resultados:

```bash
git rev-parse HEAD
git status --short
```

## Único cálculo largo vigente: convergencia temporal de una celda

Ejecutar manualmente después de revisar el plan. Compara SSPRK3 con 160, 320 y
640 pasos contra una referencia del mismo método con 1280 pasos. El propósito
es verificar si el error de energía disminuye como requiere el criterio temporal;
este ensayo no cierra por sí solo la etapa 2.

Plan: `docs/implementation/stage2/review_20260922/manual_time_plan.json`.
Duración estimada: unos 26 minutos, un hilo de CPU; reservar 2 GiB de RAM
(reserva estimada; no se ha medido el pico de memoria de este lote).
La estimación puede variar. El lote conserva cada resultado y se detiene ante
el primer fallo. No reintenta ni sobrescribe una carpeta de resultados existente.

Vista previa opcional, sin integrar las ecuaciones:

```bash
"$NOTEBOOK_PY" sandbox/stage2_cells/review_20260922/run_manual_time_batch.py --plan docs/implementation/stage2/review_20260922/manual_time_plan.json --output-root tmp/stage2_ssp_time_20260922 --dry-run
```

Comando para el cálculo completo, en la misma terminal preparada arriba:

```bash
"$NOTEBOOK_PY" -u sandbox/stage2_cells/review_20260922/run_manual_time_batch.py --plan docs/implementation/stage2/review_20260922/manual_time_plan.json --output-root tmp/stage2_ssp_time_20260922 --execute
```

Las salidas serán `one_ssp_160.json`, `one_ssp_320.json`, `one_ssp_640.json`,
`one_ssp_1280.json` y `one_ssp_time_assessment.json`, más NPZ, registros y recibos.
Conservar la carpeta `tmp/stage2_ssp_time_20260922/` completa e indicar «reinicia»
cuando termine o se detenga. El agente revisará sus resultados antes de proponer
otra ejecución. No repetir el comando sobre la misma carpeta.

Los cálculos previstos de más de cinco minutos los ejecuta el usuario. No se
divide ni reintenta un trabajo para eludir ese límite. Si una prueba limitada
vence su tiempo, su resultado queda incompleto y se revisa antes de continuar.

No reutilizar los lanzamientos del histórico ni las instrucciones de ejecución
de manifiestos anteriores como una continuación automática de esta revisión.
"""


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "docs/GEMINGA_COMMANDS.md")
    parser.add_argument("--draft", type=Path, default=DEFAULT_DRAFT)
    parser.add_argument("--expected-source-sha256", default=EXPECTED_SOURCE_SHA256)
    args = parser.parse_args()
    source = args.source.resolve()
    draft = args.draft.resolve()
    archive = ARCHIVE.resolve()
    active = (ROOT / "docs/GEMINGA_COMMANDS.md").resolve()
    if draft in (source, archive, active) or archive == source:
        raise SystemExit("Refusing to overwrite the source, active mirror, or archive with a draft.")
    if not draft.is_relative_to(ROOT) or not archive.is_relative_to(ROOT):
        raise SystemExit("Outputs must remain inside the repository.")
    original = source.read_bytes()
    original_hash = sha256(original)
    if original_hash != args.expected_source_sha256.lower():
        raise SystemExit(f"Source changed: observed SHA256 {original_hash}; review before preparing a draft.")
    if archive.exists() and archive.read_bytes() != original:
        raise SystemExit("Archive already exists with different bytes; preserving it and stopping.")
    archive.parent.mkdir(parents=True, exist_ok=True)
    if not archive.exists():
        archive.write_bytes(original)
    prepared = DRAFT.encode("utf-8")
    if "screen" in DRAFT.lower():
        raise SystemExit("Unexpected terminal-session instruction in active draft.")
    draft.parent.mkdir(parents=True, exist_ok=True)
    draft.write_bytes(prepared)
    if source.read_bytes() != original or archive.read_bytes() != original:
        raise SystemExit("Source/archive integrity verification failed.")
    report = {
        "schema": "pysnspd.command_notebook_cleanup_draft.v1",
        "status": "DRAFT_ONLY_NO_ACTIVE_OR_REMOTE_WRITE",
        "source": str(source),
        "source_sha256": original_hash,
        "archive": str(archive),
        "archive_sha256": sha256(archive.read_bytes()),
        "draft": str(draft),
        "draft_sha256": sha256(prepared),
        "original_bytes": len(original),
        "draft_bytes": len(prepared),
        "original_lines": len(original.splitlines()),
        "draft_lines": len(prepared.splitlines()),
        "source_unchanged": True,
        "original_terminal_session_mentions": original.lower().count(b"screen"),
        "draft_terminal_session_mentions": prepared.lower().count(b"screen"),
        "active_long_computation_commands": 1,
    }
    report_path = draft.with_name("GEMINGA_COMMANDS_cleanup_provenance.json")
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
