# Geminga: comandos vigentes

Actualizado: 22 de septiembre de 2026. Cuenta `jdiaz`, sin administrador.
Repositorio `/home/jdiaz/pysnspd`; base Git `031991e` más revisión de trabajo.

## Resultado recibido y decisión

El lote `tmp/stage2_ssp_two_time_20260922/` terminó en 48,54 minutos:
**98/98 controles PASS**. La convergencia temporal de una y dos celdas pasa
para el integrador SSP anterior. Esas corridas quedan conservadas.

Una prueba adicional de transporte aislado encontró una ocupación negativa
de -2,5e-323 causada por redondeo de inventarios subnormales. La nueva protección
aritmética corrige ese caso sin clipping ni reparación energética y conserva
bit a bit el cálculo ordinario ensayado. Sus pruebas ligeras se han revisado;
falta acreditar su propia convergencia y las mallas completas.

**El antiguo plan de sólo mallas `remaining_mesh_plan.json` está bloqueado.**
No ejecutarlo. El nuevo plan vuelve a medir cada gate dinámico con la versión
protegida. Esto satisface la validación independiente de cada versión; los
resultados anteriores no se usan para aprobar el nuevo integrador.

Histórico íntegro de instrucciones anteriores:
[GEMINGA_COMMANDS_before_mesh.md](/home/jdiaz/pysnspd/docs/implementation/stage2/closure_prep_20260922/GEMINGA_COMMANDS_before_mesh.md).
Es evidencia histórica, no una lista para volver a ejecutar.

## Preparación de la terminal

```bash
cd /home/jdiaz/pysnspd
NOTEBOOK_PY=/home/jdiaz/.conda/envs/snspd/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Cálculo largo pendiente: validación completa del integrador protegido

Propósito: comprobar soporte de poblaciones, balances, precisión temporal y
convergencia en las dos mallas energéticas antes de cerrar la etapa 2.
**35 tareas, con 19 trayectorias nuevas**, en este orden:

1. Equilibrio, vacío fonónico y electrones dispersos: trayectorias nuevas,
   invariantes y campos/soporte; incluye el equilibrio con Gamma finito.
2. Una celda y dos celdas: 160/320/640 pasos frente a su referencia propia de
   1280. La pareja TWO320/640 además debe pasar el presupuesto de 2,5e-5.
3. Fonones 2049/4097 y electrones 1260/2520, cada malla con 320/640 pasos.
   Se comprueba el error temporal de cada pareja antes de comparar tres mallas
   a 320 pasos: 1025/2049/4097 y 630/1260/2520.
4. Campos y soporte sobre las trayectorias finalmente seleccionadas.

Los límites registrados no cambian. Entradas Debye sintéticas; no se evalúan
tasas absolutas de NbN, circuitos ni un transiente espacial del detector.

**Estimación: 9 a 15 horas, un hilo; reservar 12 GiB de RAM.** El pico de RAM y
el coste de la aritmética extendida sobre la malla de 2520 estados no están
medidos; el intervalo es una estimación, no una garantía. Si se activa mucho
esa rama, puede tardar más. No hay procesos de fondo ni esperas del agente.

Vista previa opcional, sin dinámica y sin crear resultados:

```bash
"$NOTEBOOK_PY" sandbox/stage2_cells/closure_prep_20260922/run_guarded_acceptance_batch.py --plan docs/implementation/stage2/closure_prep_20260922/guarded_acceptance_plan.json --output-root tmp/stage2_guarded_acceptance_20260922 --dry-run
```

**Único comando largo vigente**, en la misma terminal:

```bash
"$NOTEBOOK_PY" -u sandbox/stage2_cells/closure_prep_20260922/run_guarded_acceptance_batch.py --plan docs/implementation/stage2/closure_prep_20260922/guarded_acceptance_plan.json --output-root tmp/stage2_guarded_acceptance_20260922 --execute
```

Salidas: `tmp/stage2_guarded_acceptance_20260922/`, con 19 trayectorias JSON/NPZ,
estados iniciales, evaluaciones temporales/pareadas/de malla, campos/soporte,
manifiesto, recibos y logs. La carpeta no existe al entregar este plan.
Cada tarea se verifica antes de continuar. Al primer fallo se detiene y
conserva todo: no modifica parámetros, reintenta ni sobrescribe evidencia.
No relanzar automáticamente si se detiene. Un lote completo aún requiere
revisión conjunta de resultados antes de admitir la etapa 2 y hacer el push.

Al terminar o fallar, basta indicar la carpeta y escribir «reinicia».
No hace falta pegar el output extenso. La etapa 3 permanece preparada y sin
iniciar. No se ejecutará trabajo largo desde el agente ni se dividirá para
eludir el límite de cinco minutos.

## Comandos ligeros útiles

Verificación de guardas del lote, sin RHS (segundos):

```bash
timeout --signal=TERM --kill-after=5s 240s "$NOTEBOOK_PY" sandbox/stage2_cells/closure_prep_20260922/run_guarded_acceptance_batch.py --self-test
```

Auditar otra vez la recepción TWO, sin reintegrar (segundos):

```bash
timeout --signal=TERM --kill-after=5s 240s "$NOTEBOOK_PY" sandbox/stage2_cells/closure_prep_20260922/audit_two_cell.py --output tmp/two_cell_audit_repeat.json
```

Verificar la entrega exacta por hashes (segundos):

```bash
"$NOTEBOOK_PY" sandbox/stage2_cells/verify_delivery.py
```

No ejecutar diagnósticos sobre sus archivos congelados. Evidencia y alcance:
`docs/implementation/stage2/closure_prep_20260922/closure_status.json`.
