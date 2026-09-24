# Referencia térmica espacial: interfaz e integración acotada

Se puede incorporar una energía espectral espacial sin sustituir ahora el
integrador ni la malla de producción. La implementación inicial está en
`pysnspd/experimental/thermal_spatial_usadel.py`: opera sobre áreas nodales,
aristas orientadas y conductancias geométricas. El constructor cartesiano
permite comprobarla antes de adaptar la geometría Delaunay–Voronoi existente.
«Exacta» aquí significa derivadas de una misma energía discreta con el mismo
corte de Matsubara; no ausencia de errores de malla o truncamiento.

## Energía y variables espectrales

Se usan d=Δ/(kBTc), t=T/Tc, εₙ=2πt(n+1/2), y coordenadas X=x/ℓ₀, con
ℓ₀²=ℏD/(2kBTc). En cada nodo, una variable compleja u parametriza

```math
f=\frac{u}{\sqrt{1+|u|^2}},\qquad
g=\frac1{\sqrt{1+|u|^2}},\qquad
\mathbf n=(\Re f,\Im f,g).
```

Así se conserva |n|=1 y g>0 en frecuencias positivas; u es regular cuando
d=0. No se emplea qδ ni una división por la amplitud del condensado.

Sean mᵢ las áreas adimensionales, cᵢⱼ≥0 las conductancias geométricas y
Uᵢⱼ=exp(−iαᵢⱼ). La energía cruda con N frecuencias es

```math
F_N=\sum_i m_i|d_i|^2\ln t+
2\pi t\sum_{n=0}^{N-1}\left\{
\sum_i m_i\left[\frac{|d_i|^2}{\epsilon_n}
+2\epsilon_n(1-g_{ni})-2\Re(d_i^*f_{ni})\right]
+\sum_{(i,j)}c_{ij}\left[|f_{ni}-U_{ij}f_{nj}|^2+(g_{ni}-g_{nj})^2\right]
\right\}.
```

Su variación espectral reproduce la ecuación radial ya ensayada cuando se
impone esa simetría. No se añade un laplaciano K₀ externo: la respuesta espacial
ya está incluida. Para una película uniforme, la escala de esta energía 2D es
N₀(kBTc)²·espesor·ℓ₀².

En la solución espectral estacionaria, el teorema de la envolvente da el
gradiente cartesiano complejo integrado

```math
G_i=2m_i d_i\ln t+4\pi t m_i\sum_n(d_i/\epsilon_n-f_{ni}).
```

La derivada respecto de α es −4πt cᵢⱼΣₙIm(f*ₙᵢUᵢⱼfₙⱼ). Se devuelve por
separado `link_derivative_alpha` y `current_bar=-link_derivative_alpha`,
orientada de cola a cabeza. Esta distinción evita mezclarla con la convención
histórica que usa el exponente de enlace con signo contrario. Para pasar la
corriente a SI se multiplica `current_bar` por la escala energética y 2e/ℏ.

## Resolución e interfaz disponibles

`ThermalGraph(area_weights, edges, conductance, ...)` acepta conductancias
nulas como aristas inactivas. `rectangular_graph(nx,ny,L,W)` usa áreas
trapezoidales y medidas transversales reducidas a la mitad en el borde. No
introduce diagonales artificiales ni condiciones periódicas.

`solve_frequency(graph,d,epsilon,alpha,...)` resuelve una frecuencia independiente.
El residuo es Rᵢ=Hzᵢuᵢ−Hxyᵢ, con
Hzᵢ=mᵢε+Σⱼcᵢⱼgⱼ y Hxyᵢ=mᵢdᵢ+ΣⱼcᵢⱼUᵢⱼgⱼuⱼ. Su Jacobiano real contiene
HzᵢI en la diagonal y, fuera de ella,

```math
c_{ij}\left[-g_j\mathcal R_{ij}
+g_j^3(\mathcal R_{ij}u_j-u_i)u_j^T\right],
```

donde Rᵢⱼ rota las dos componentes de u. Newton usa el bloque libre disperso
y una búsqueda Armijo que reduce la energía. Se registran residuo, iteraciones
y reducciones. Un paso sin descenso, una búsqueda fallida o falta de convergencia
detienen ese modo: no se cambia de algoritmo ni se recorta su solución.

`spectral_energy_gradient` y `spectral_residual_jacobian` exponen las derivadas
para comprobaciones independientes. `evaluate_thermal` agrega soluciones de
frecuencias consecutivas y devuelve energía, fuerza y corriente de la misma
suma finita. Comprueba que las soluciones pertenecen al grafo, campo y enlaces
de la evaluación; no reutiliza silenciosamente un espectro de otro estado.

## Bordes y reutilización del repositorio

Los nodos espectrales libres tienen el borde natural sin flujo de esta energía.
Los valores espectrales impuestos se pasan como `fixed_nodes` y `fixed_u`.
En una variación de d se conservan esos valores, o se incluye explícitamente su
dependencia y trabajo de reacción. La identidad de calibre debe contar también
la reacción de fase de las variables espectrales fijadas. Esto todavía no es
un reservorio de corriente ni el contacto completo del detector.

Los adaptadores futuros pueden tomar `node_area_m2`, `edge_length_m` y
`dual_face_length_m` de `pysnspd/mesh/operators.py`: m=área/ℓ₀² y c=longitud
dual/longitud primal. La incidencia cola(+)/cabeza(−), los tags de borde de
`mesh/edges.py` y el álgebra de `experimental/electrical_ports.py` ya son útiles.
La corriente integrada puede entrar en ese puerto sin reconstruirla mediante
una interpolación independiente.

Hay un detalle que impide considerar trivial ese adaptador: el constructor
actual de duales reemplaza longitudes nulas o inválidas por media arista primal.
Una diagonal de una triangulación cartesiana puede tener dual exactamente
nulo. Se debe auditar la procedencia geométrica y conservar ese cero cuando
corresponda; no declarar una referencia espacial por haber importado pesos
positivos. El adaptador de producción no se incluye en esta primera entrega.

## Suma espectral, coste y alcance

La implementación inicial devuelve sólo F_N, sin una cola escondida. Si se
añade una cola espacial, debe añadirse primero su energía y derivarse tanto
fuerza como corriente. Por ejemplo, el primer término es
2πt S₂Σc|dᵢ−Uᵢⱼdⱼ|², con S₂=Σₙ≥Nεₙ⁻²; no basta agregar una corrección a la
fuerza radial. Su validez requiere frecuencias altas frente a las escalas
espaciales realmente presentes. La comparación de cortes permanece separada
de la comparación física entre cierres.

Un primer piloto 2D de 9×5 nodos, 76 aristas, dos perfiles y ocho frecuencias
supone 16 problemas no lineales de 90 incógnitas reales. Cada evaluación es
O(nodos+aristas) y cada modo tiene su sistema disperso. Las siete pruebas
actuales, sobre 20 nodos, tomaron 0,201 s de cálculo local; esto no es todavía
un benchmark del piloto físico. Una expectativa inicial de segundos a decenas
de segundos debe verificarse con un límite de 240 s, sin reinicios automáticos.

Los pares perfil/frecuencia pueden compartir el pool de
`sandbox/stage4_core/parallel_runtime.py`, respetando presupuesto efectivo CPU,
memoria, afinidad y BLAS de un hilo. No se necesita paralelismo anidado. Si un
caso exige más de cinco minutos, se entrega el comando al usuario conforme a
la política vigente.

Este módulo es una referencia **térmica y estática**. No proporciona aún una
DOS retardada, transporte no térmico, movilidad KWT calibrada o una evolución
autoconsistente de Δ. No se inserta directamente como fuerza en la cinética de
etapa 2 sin revisar la compatibilidad energética correspondiente, ni modifica
el circuito de la memoria o el solver de producción.
