# Malla dual preparada para el siguiente control suave

**Usar `resampled/mesh.npz` para el siguiente ensayo.** Una sola corrección del
muestreo de borde eliminó los avisos y el exceso de área de la primera malla.
La primera preparación descrita más abajo se conserva íntegra como evidencia.

| Medida geométrica | Primera preparación | Muestreo de borde corregido |
|---|---:|---:|
| Nodos | 2126 | 1712 |
| Aristas | 6187 | 4941 |
| Nodos por contacto | 32 | 33 |
| Avisos del constructor | 13 | 0 |
| Circuncentros fuera del rectángulo | 8 | 0 |
| Exceso de área dual | 0,2345 % | 0 al redondeo medido |
| Diferencia relativa con laplaciano heredado | 2,04×10⁻¹⁴ | 1,60×10⁻¹⁴ |

El ajuste usa `boundary_points=196`: con la asignación publicada de
[pyTDGL `box`](https://py-tdgl.readthedocs.io/en/stable/_modules/tdgl/geometry.html)
resultan 64 intervalos longitudinales y 32 transversales, todos de 2,5 nm,
la mitad de la longitud de arista objetivo del interior. La forma y las
dimensiones físicas no cambian. No se probaron semillas aleatorias ni una
familia de tolerancias. Esta segunda preparación tardó 0,114 s.

La primera triangulación tenía ocho circuncentros fuera del rectángulo,
hasta 1,12 nm más allá del borde; su suma de áreas triangulares sí coincidía
con el rectángulo. La corrección elimina esos circuncentros exteriores y el
exceso de área dual. Este diagnóstico geométrico se obtuvo directamente de
los vértices de los triángulos; consta en `geometry_comparison.json`. No se
modificó el algoritmo de áreas ni se recortaron sus resultados.

## Primera preparación conservada

La malla fue generada por la ruta `Device.make_mesh` y `Mesh.from_triangulation`
ya presente en producción. No se resolvieron espectros ni transientes. El
archivo `mesh.npz` conserva exactamente sus áreas y caras duales; el adaptador
no recalcula, mezcla ni recorta esas medidas.

El rectángulo tiene ancho **80 nm** y largo **160 nm**, con arista máxima de
5 nm. La relación largo/ancho 2 es provisional para este control. No es una
conclusión sobre el dominio espacial adecuado del experimento de Korzh.

- 2126 nodos, 4062 triángulos y 6187 aristas.
- 32 nodos de contacto en cada extremo longitudinal; los laterales permanecen
  como bordes naturales del grafo.
- ℓ0 = 4,698490895830414 nm, con D = 0,5 cm²/s y Tc = 8,65 K.
- Coordenadas SI: x entre 0 y 160 nm, y entre −40 y 40 nm.
- `coordinates_bar` traslada el origen al centro (80 nm, 0) y divide por ℓ0.
- `smooth_profile` es cos²[π(x−L/2)/L] cos²(πy/W), adimensional entre 0 y 1.
  Es la forma de una perturbación suave de prueba; no tiene energía fotónica,
  ancho de cascada ni temperatura inicial asignados.

La preparación tardó **0,160 s**. Tres pruebas del adaptador pasaron en 0,011 s.
La aplicación del laplaciano al perfil preparado coincide con la rutina
heredada `build_laplacian` hasta una diferencia relativa de **2,04×10⁻¹⁴**.
Las pruebas también comprueban el uso de las mismas medidas por el operador
FV y la invariancia al reescalar conjuntamente longitudes y ℓ0.

El constructor heredado emitió 13 avisos de ordenación de celdas de borde,
conservados en `preparation.log`. La suma de áreas es
1,2830018413×10⁻¹⁴ m², un **0,2345 %** mayor que el área geométrica nominal.
Todas las áreas y conductancias resultaron positivas. Esta diferencia queda
registrada para el presupuesto espacial y los balances; no se ajustó el área
ni se cambió producción para hacerla desaparecer. La identidad con el operador
heredado prueba la transferencia del mismo operador, no exactitud geométrica
universal ni convergencia del futuro transiente.

El recibo registra hashes, parámetros y coste. `execution_resources.json`
contiene la topología, afinidad, cuota y memoria consultadas en Geminga:
la preparación usó un solo proceso fijado a un CPU, con bibliotecas numéricas
limitadas a un hilo. El campo `seed=12345` se conserva en los parámetros, pero
esta ruta concreta de MeshPy no consume explícitamente una semilla aleatoria.

Para un gráfico de esta preparación, el rótulo correcto es **«Forma espacial
adimensional de la perturbación de prueba»**, con x−L/2 e y en nm y barra de
color «perfil, sin unidades». Si se muestra la malla, identificar aristas
Delaunay, caras Voronoi y contactos; no llamarla mapa de temperatura o hotbelt.

Comando reproducible de la **primera preparación**, con una carpeta nueva:

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/prepare_dual_mesh.py \
  --output-root tmp/stage4_dual_mesh_manual
```

El lanzador de campaña que utilice esta malla debe aplicar su presupuesto
común de afinidad, memoria y trabajadores antes de resolver los espectros.
