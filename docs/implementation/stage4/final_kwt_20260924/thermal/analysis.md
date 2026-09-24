# Ensayo térmico final sobre malla dual

**La trayectoria completó 1 ps en 153.81 s y cumplió todos los márgenes registrados.** Se admitieron 40/40 comparaciones temporales a tiempo positivo, además de ocho comparaciones coincidentes del estado inicial. Se usó el paso Euler/KWT real de la memoria, de primer orden.

La malla tiene 1712 nodos sobre un rectángulo de 160 × 80 nm. Se estudiaron perturbaciones suaves de amplitud y fase, con temperatura fija de 0,9 K, contactos superconductores de equilibrio en los extremos y paredes laterales aislantes. Son sondas del núcleo térmico; no representan depósito fotónico, hotbelt ni evolución cinética de poblaciones.

| Sonda | Magnitud | Mayor diferencia temporal (% de la señal inicial) | Tiempo (ps) |
|---|---|---:|---:|
| Amplitud | Campo del condensado | 0.00150509 | 1 |
| Amplitud | Corriente espectral | 0 | 0 |
| Amplitud | Fuerza por área | 0.00151068 | 1 |
| Amplitud | Fuerza de fase por área | 0 | 0 |
| Fase | Campo del condensado | 0.148423 | 0.1 |
| Fase | Corriente espectral | 0.20239 | 0.1 |
| Fase | Fuerza por área | 0.0150235 | 0.1 |
| Fase | Fuerza de fase por área | 0.211549 | 0.1 |

El denominador de cada comparación temporal es la mayor norma inicial de **esa misma magnitud** entre las dos sondas. El margen registrado es 2% más un piso absoluto de 10⁻⁷ en la norma correspondiente; no se divide por una cola tardía extinguida.

| Trayectoria | Pasos aceptados | Máximo residuo del balance (% del exceso inicial) |
|---|---:|---:|
| primary_amplitude | 581 | 0.00605432 |
| primary_angular_phase | 581 | 1.14673 |
| refined_amplitude | 1162 | 0.00302598 |
| refined_angular_phase | 1162 | 0.577085 |

El balance integra la disipación KWT y Joule normal en todos los pasos aceptados, con tiempo físico en ps. Se compara con la caída de la energía libre respecto al equilibrio uniforme. No se observó aumento de energía entre pasos. Este balance térmico no constituye una prueba de conservación de energía interna del sistema de poblaciones.

Se verificaron 893,440 predicciones espectrales mediante el residuo no lineal exacto; ninguna necesitó la corrección Newton disponible. El máximo residuo fue 3.961e-08, inferior al límite registrado de 1e-07. La predicción reutiliza factorizaciones, sin sustituir el modelo por su linealización.

Se emplearon 27 trabajadores y un coordinador con un hilo numérico por proceso, dejando dos núcleos físicos completos libres. La ejecución completa está guardada; no queda pendiente repetir este ensayo ni ejecutar el comando histórico del piloto.

## Figuras y significado de las magnitudes

- `figures/01_evolucion_y_refinamiento.png`: norma del campo complejo respecto al equilibrio, normalizada por su propia sonda inicial; diferencias entre pasos temporales con el denominador registrado.
- `figures/02_balance_integrado.png`: energía libre excedente, disipación acumulada y residuo integrado durante toda la trayectoria, cada una normalizada por el exceso energético inicial de su sonda.
- `figures/03_campos_espaciales.png`: módulo del cambio del condensado sobre la malla, en porcentaje del gap uniforme; coordenadas físicas en nm y paneles separados del error temporal.

Aquí d = Δ/(kBTc), mᵢ = Aᵢ/ell0² y ||z||ₘ² = Σmᵢ|zᵢ|². Para corrientes de arista, ||I||² = Σ|Iᵢⱼ|²/cᵢⱼ. La fuerza de fase es Im(conj(d)G)/m, una fuerza variacional adimensional, no un torque mecánico.

La reproducción leyó los datos extraídos y verificó sus hashes; no ejecutó nuevas soluciones físicas. Los valores de `refinement.json` se recuperaron exactamente a partir de los campos guardados.

```bash
python -m sandbox.stage4_core.analyze_dual_kwt_results --folder docs/implementation/stage4/final_kwt_20260924/thermal
```
