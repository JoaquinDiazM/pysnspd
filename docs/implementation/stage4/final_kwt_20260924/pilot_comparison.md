# Comparación de los pilotos KWT sobre malla dual

**Actualización posterior al piloto:** la trayectoria completa de 1 ps terminó en 153.81 s y cumplió sus criterios registrados. La estimación y el consejo de entregar un comando largo se conservan aquí como contexto histórico; ya no queda una ejecución térmica pendiente. Los resultados finales están en `thermal/analysis.md`.

Los dos pilotos avanzaron hasta **0,001 ps** en la misma malla de 1712 nodos, con el mismo paso KWT heredado y 256 frecuencias. Se mantuvieron temperatura, contactos, estado inicial y tolerancia espectral 10⁻⁷. El predictor sólo cambia el punto de partida de la consulta espectral: cada resultado pasa por el residuo no lineal exacto y, si no cumple, se aplica Newton.

| Medida | Newton sin predictor | Predictor verificado |
|---|---:|---:|
| Tiempo total | 29.571 s | 4.661 s |
| Preparación y evaluación inicial | 13.737 s | 4.404 s |
| Avance aceptado y registro | 15.834 s | 0.258 s |
| Consultas que necesitaron Newton | 2560 | 0 |
| Predicciones con residuo admitido | 0 | 2560 |
| Máximo residuo espectral | 3.229e-10 | 3.961e-08 |

La mayor diferencia entre campos finales es **1.481e-06% de la perturbación inicial de su propia sonda**. Se compara d=Delta/(kB Tc), con norma sqrt(sum m|delta d|²), donde m es el área dual en unidades de ell0². Las dos resoluciones temporales usan los mismos pasos entre ambas implementaciones.

El máximo residuo del balance integrado del piloto optimizado representa **0.02077% de su exceso inicial de energía libre**. No se detectó aumento de energía entre pasos aceptados. Esto verifica el piloto de 0,001 ps; no certifica todavía el horizonte de 1 ps ni el balance de energía interna de poblaciones.

Se emplearon 27 trabajadores y un coordinador, todos con un hilo numérico, sobre 28 CPU lógicas de 32. Quedaron dos núcleos físicos completos libres. La reserva de memoria de 30 GiB estuvo por debajo del 90% disponible. El inventario y la afinidad se vuelven a comprobar en cada ejecución.

La aceleración total del piloto fue 6.34×; la del bloque de avance y registro fue 61.48×. La extrapolación inicial fue de 4.3 min y se trató como una estimación provisional basada en dos subpasos, pues podían aparecer nuevas correcciones Newton.

Las fuentes originales y su prueba permanecen idénticas a los archivos archivados del primer piloto. `pilot_comparison.json` contiene las diferencias por sonda, normas, residuos, inventarios y hashes de entrada. El análisis no ejecutó nuevas soluciones físicas.

Reproducción:

```bash
python -m sandbox.stage4_core.dual_kwt_compare_pilots --folder docs/implementation/stage4/final_kwt_20260924
```
