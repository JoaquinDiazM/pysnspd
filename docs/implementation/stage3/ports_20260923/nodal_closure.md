# Cierre de desarrollo de 3A: ensayos nodales registrados

La campaña nodal permite cerrar 3A en el alcance estático registrado y continuar el desarrollo de bordes. No cierra la etapa 3 completa ni certifica precisión espacial universal. Los resultados y las tolerancias originales se conservan.

La auditoría de datos guardados dio **PASS_SAVED_DATA_AUDIT**: ocho fuentes, 18 casos, 76 archivos y 336 nodos espaciales comprobados. El lote registró 7.845 minutos. No se evaluaron espectros, fuerzas ni trayectorias nuevas.

Las diferencias finitas de fuerza y corriente se ejecutaron en los seis casos de 16 celdas y en el piloto reutilizado de fase térmica de ocho celdas. Sus errores absolutos máximos son 1.347e-11 y 8.537e-12. Los otros 11 resultados conservan expresamente que no repitieron diferencias finitas. La covariancia de fase y el oráculo uniforme discreto alcanzan diferencias máximas de 1.794e-15 y 6.029e-16.

Las 11 respuestas no nulas tienen diez estimaciones relativas aceptadas y una sin certificado relativo. La mayor estimación es 0.613214%, frente al objetivo registrado de 1%. Corresponde a la corriente inducida por la modulación térmica de amplitud; el orden observado es 0.608. El estimador conservador incorpora esa razón medida, pero sigue siendo una extrapolación condicional; no demuestra orden cuatro de la corriente.

| Caso | Energía de exceso | Fuerza inducida | Corriente inducida |
|---|---:|---:|---:|
| weak_amplitude_thermal | 0.000533496% | 0.000573521% | 0.613214% |
| weak_phase_thermal | 0.0426285% | 0.0121958% | 0.0411617% |
| weak_phase_nonthermal | 0.0423539% | 0.0122075% | 0.0404925% |
| weak_amplitude_vacuum | Bajo piso declarado; sin pase relativo | 3.12442e-05% | Cero analítico; pase absoluto |

En amplitud/vacío, la diferencia de energía entre 16 y 32 celdas es 1.835e-9, inferior al piso declarado de 1e-8. Se conserva `ROUND_OFF_LIMITED`; no significa que se haya demostrado el límite de redondeo de la máquina ni se convierte en un pase relativo. La corriente exactamente nula pasa sólo su control absoluto. Los uniformes se cotejan con su oráculo discreto, sin cocientes 0/0.

El error del gradiente frente a la referencia analítica disminuye en 8, 16 y 32 celdas: 0.400155%, 0.0260627%, 0.00164582%. El valor fino cumple el objetivo de 1%.

El menor valor propio local D.36 es 1.57079633, frente a una incertidumbre de diferenciación máxima de 2.793e-05. Los 12 puntos representativos comparados entre 630 y 1260 nodos ocupacionales conservan el signo positivo. La variación máxima de matriz es 4.054e-09; la variación escalada máxima de momentos es 5.394e-09. El control negativo conserva el valor propio -0.30661962 y su rechazo esperado. Estas comprobaciones no prueban estabilidad dinámica, todo el Hessiano nodal ni el espectro continuo.

Se reutiliza la evidencia periódica íntegra. Lo siguiente requiere controles propios de bordes abiertos, intercambio con reservorios, continuidad de corriente y balance de trabajo/potencia del circuito. Después corresponde el ensayo espacial en el tiempo y su refinamiento sobre observables del dispositivo. No hace falta repetir esta campaña sin cambios para iniciar ese trabajo.

Fuentes y trazabilidad: registro nodal `../nodal_20260923/registration.json`; resultados íntegros `raw/nodal_campaign/`; inventario SHA256, dictámenes y reproducción numérica en `nodal_audit.json`; script `sandbox/stage3_spatial/ports_20260923/audit_nodal.py`.
