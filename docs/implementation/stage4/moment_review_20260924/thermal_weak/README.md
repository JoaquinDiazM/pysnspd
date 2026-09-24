# Operador térmico débil: preparación completada

La campaña ya terminó y sus resultados permiten pasar a la evolución temporal
afín de todos los nodos. **No es necesario repetir esta campaña.** Se conserva
el núcleo radial aceptado, su contorno espectral R12, la temperatura de 0,9 K
y los coeficientes KWT heredados. Producción permanece intacta.

Se resolvieron 256 frecuencias de Matsubara y 2304 raíces espectrales en
**32,94 s**, usando 27 trabajadores y un coordinador, un hilo numérico por
proceso y dos núcleos físicos libres. El piloto de coste duró 4,35 s.

| Comprobación | Resultado |
|---|---:|
| Error relativo del Hessiano, paso 0,01 | ≤7,08×10⁻⁷ |
| Error relativo de la corriente, paso 0,01 | ≤2,32×10⁻⁶ |
| Error relativo del RHS KWT y potencial, paso 0,01 | ≤4,86×10⁻⁵ |
| Continuidad nodal | ≤1,67×10⁻¹⁶ |
| Residuo del balance de energía libre | ≤1,26×10⁻¹² |
| Identidad gauge linealizada | ≤1,64×10⁻¹⁰ |
| Checkpoints con SHA verificado | 256/256 |

Los errores de las diferencias finitas disminuyen aproximadamente por cuatro
al dividir el paso por dos. La derivada de movilidad y el potencial de la
referencia aportan 1,10–1,13 % del RHS ensayado; ambos términos están incluidos.
Los dos perfiles de prueba verifican aplicaciones del operador: no reducen
el sistema físico a dos variables.

La validación corresponde a la energía libre térmica, con temperatura y
contornos fijos. No acredita una trayectoria no lineal, el balance de energía
no térmica, un transiente fotónico ni una señal de detector. La etapa 4 sigue
abierta a esos acoplamientos.

## Comprobantes y reutilización

- [Resumen numérico](raw/summary.json), [identidad y recursos](raw/identity.json)
  y [progreso de la ejecución](raw/progress.jsonl).
- [Campos y aplicaciones del operador](raw/full_node_operator_checks.npz).
- [Verificación de integridad](verification_receipt.json) y
  [piloto de coste](pilot_summary.json).
- [Registro de las cinco pruebas unitarias](unit_tests.log) y su
  [recibo con hashes](unit_tests_receipt.json): todas pasan, en 0,462 s.
- [Plan ejecutado](plan.json) y [prueba de que el dry-run no calcula](dry_run_guard_proof.json).

Las 256 matrices espectrales, soluciones y metadatos completos ocupan
370 MB y permanecen en Geminga, en
`/home/jdiaz/scratch/stage4_thermal_weak_20260924`. El driver temporal reutiliza
esas matrices para construir una vez sus factores LU y aplicar el Jacobiano
en todos los nodos; no necesita repetir las 2304 raíces. Debe conservar la
deriva de referencia y distinguir el balance cuadrático de la evolución afín
del balance no lineal de energía libre.
