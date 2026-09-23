# Revisión del postproceso con reservorios

El postproceso queda verificado como ensayo instantáneo: 12 fuentes, 14 artefactos de entrada y los cuatro artefactos de salida de los casos coinciden con sus hashes. Los campos, poblaciones y masas usadas permanecen idénticos a los archivados. No se consultó un espectro nuevo ni se integró un paso temporal.

La recuperación E=theta[log(1-p)-log(p)] es la inversión de la preparación de Fermi–Dirac registrada. Se exige 0<p<0.5, energías positivas y ordenadas, y reproducción del momento de calor original. La revisión recompone además el depósito nuevo y CM9 directamente de los arrays guardados. Esta recuperación no permite usar esas energías para otros campos o poblaciones fuera del estado archivado.

La carga de fase procede de Iref, no del gradiente discreto. La reacción radial impone radio fijo dentro del solve KWT y mantiene explícito su trabajo, numéricamente nulo. El potencial se resuelve con Iref y el circuito usa Ib=Is=Iref, vc=0. El trabajo material de reservorio entra una vez en CM9 y no se suma al calor electrónico.

| Magnitud | Hélice | Perturbación suave |
|---|---:|---:|
| QDelta con extremos libres | 1.01145943 | 1.01153896 |
| QDelta con carga de reservorio | 6.19215888e-07 | 8.01522847e-05 |
| QDelta terminal con carga | 1.77720022e-07 | 1.77720022e-07 |
| Velocidad material máxima | 0.0011960404 | 0.0140151015 |
| Vdev [microvoltios] | 0.904496123 | 0.915340354 |

QDelta y velocidad usan las unidades normalizadas del registro. El trabajo de reservorio es -6.82937323 pW en ambos estados. CM9 recompuesta deja residuos de orden1e-27 W. El radio terminal cambia menos de1.4e-18 por unidad de tiempo normalizada. El defecto normal en la primera cara interior es aproximadamente -0.04193495% de Iref: se conserva como diagnóstico, no se fuerza a cero.

El descenso de QDelta revela la importancia de imponer una condición terminal compatible. El circuito también cambió de Is=0.99Iref a Is=Iref; por ello la comparación de voltaje o calor Joule no es una medida aislada de mejora numérica.

Recomiendo las seis instantáneas cargadas para estudiar sensibilidad espacial de QDelta, de la respuesta interior inducida y del defecto normal adyacente. El control libre permanece archivado. No se obtiene de ello un transiente, radio dinámico a_b(Is), intercambio cinético de reservorio ni admisión completa de D.27. La corriente normal externa nula es parte del contrato impuesto; la cara interior no certifica por sí sola esa traza continua.

El primer intento, fallido antes de KWT por comparar masas Linux/Windows bit a bit, permanece íntegro. La diferencia geométrica fue2.58e-15 relativa; el nuevo postproceso valida geometría dentro de64eps y usa exactamente las masas archivadas. Los criterios físicos no cambiaron. El registro declara la exploración algebraica anterior y no se presenta como un prerregistro ciego.

Reproducción sólo por aritmética guardada: `python sandbox/stage3_spatial/coupled_20260923/audit_saved_reservoir.py`. Rechaza sobrescribir esta revisión.
