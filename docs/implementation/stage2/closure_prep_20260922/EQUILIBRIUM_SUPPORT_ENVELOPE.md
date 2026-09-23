# Envolvente analítica para la absorción fuera del soporte electrónico

Esta nota justifica una cota superior sin cuadratura para la absorción que lleva
un electrón representado a una energía superior al catálogo. Se aplica por
celda y tiempo comprobado. Las energías usan las unidades adimensionales del
ensayo; `a` es la amplitud del condensado en unidades de `Delta0`.

## Cota espectral, incluido un regulador finito

Sean `E > a`, `Gamma >= 0` y `eta > 0`. La rama retardada causal se escribe

\[
z=E+i\eta,\qquad u=z+i\Gamma c,\qquad
w=\sqrt{u^2-a^2},\qquad c=u/w,\qquad s=ia/w.
\]

Se sigue la rama con `Re u > 0`, `Im u > 0`, `Re c > 0` y `Re w > 0`.
Estas condiciones pertenecen a la selección causal, no a cualquier raíz del
polinomio espectral. En la continuación física desde `Gamma = 0`, `Re u` no
puede cruzar cero para `E > 0`: si `u = iv`, con `v > 0`, entonces
`c = v/sqrt(v²+a²)` es real positivo y la parte real de `u=z+i Gamma c`
exigiría `0=E`.

Escribiendo `u = u_R+i u_I`,

\[
\operatorname{Im}(c^2)
=\operatorname{Im}\!\left(1+\frac{a^2}{u^2-a^2}\right)
=-\frac{2a^2u_Ru_I}{|u^2-a^2|^2}\leq0.
\]

Como `Re c > 0`, esto implica `Im c <= 0`; por consiguiente
`Re u = E-Gamma Im c >= E`. La desigualdad triangular da

\[
|c|^2=\frac{|u|^2}{|u^2-a^2|}
\leq\frac{|u|^2}{|u|^2-a^2}
\leq\frac{E^2}{E^2-a^2}.
\]

La última función decrece para `E > a`. Con `N=Re c` y `R=Im s`, se tiene
`N <= |c|` y `R=a Re(1/w) >= 0`. Por tanto, para ambas energías mayores que
`E_min > a`, el factor de dispersión satisface

\[
\mathcal C_S(E,E+\Omega)
=N(E)N(E+\Omega)-R(E)R(E+\Omega)
\leq N(E)N(E+\Omega)
\leq C_{\max}:=\frac{E_{\min}^2}{E_{\min}^2-a^2}.
\]

La prueba conserva `eta` durante todos los pasos. No requiere convertir una
diferencia de regulador medida en unos puntos en una cota global. Para `a=0`
se recupera directamente `C_max=1`.

## Cota de ocupaciones sin extrapolación

Sea `E_max` la mayor energía electrónica representada y
`E_min=E_max-Omega_D > a`. Se toma `p_max` como el máximo del interpolante
electrónico lineal sobre `[E_min,E_max]`, incluyendo **ambos extremos interpolados
y todos los nudos interiores**. Se toma `n_max` sobre el soporte fonónico completo
`[Omega_lo,Omega_D]`. Los dos intervalos deben estar cubiertos por sus tablas.

Para pesos de interpolación no negativos que suman uno, la desigualdad entre
medias geométrica y aritmética asegura

\[
\prod_j p_j^{b_j}\leq\sum_j b_jp_j\leq p_{\max}.
\]

Así, el máximo del interpolante lineal domina la actividad electrónica
geométrica utilizada en los eventos. El factor de huecos del electrón externo
es a lo sumo uno para cualquier ocupación física; no se extrapola esa ocupación.
El producto Bose geométrico sin normalizar también es menor o igual a `n_max`.
Para la versión normalizada usada por el operador,

\[
n_{\rm eff}=\frac{1}{\exp\!\left[\sum_k b_k\ln(1+1/n_k)\right]-1}
\leq n_{\max},
\]

porque cada afinidad `ln(1+1/n_k)` es al menos `ln(1+1/n_max)`. Los ceros se
entienden por continuidad: un factor nulo con peso positivo impide absorción.
La tasa de absorción queda entonces dominada por `p_max*n_max`, sin recortar
ocupaciones ni introducir una población externa.

## Primitivas exactas de la envolvente

En este ensayo, `alpha2F(Omega)=lambda*(Omega/Omega_D)²` dentro del soporte.
Para cada `Omega`, sólo los electrones con
`E in [E_max-Omega,E_max]` pueden salir por absorción. Ese intervalo tiene
longitud `Omega`. Si `P` es el prefactor de reacción registrado —en las unidades
actuales, `8*pi*Delta0*t_ref/hbar`—, se obtiene

\[
\begin{aligned}
B_{\rm conteo}
&=P C_{\max}p_{\max}n_{\max}
\frac{\lambda}{\Omega_D^2}
\frac{\Omega_D^4-\Omega_{\rm lo}^4}{4},\\
B_{\rm potencia}
&=P C_{\max}p_{\max}n_{\max}
\frac{\lambda}{\Omega_D^2}
\frac{\Omega_D^5-\Omega_{\rm lo}^5}{5}.
\end{aligned}
\]

El segundo resultado incorpora la energía fonónica `Omega` de cada absorción.
Ambos son mayorantes, no estimaciones de la tasa real. Sus cocientes se comparan
con las respectivas actividades brutas registradas, manteniendo la puerta
vigente de `1e-3` y la misma convención de prefactores y conteo.

## Alcance y evidencia

El cambio relativo de la estimación numérica con órdenes 16/32 superó su
presupuesto de referencia; ese resultado fallido permanece conservado. La
envolvente aquí demostrada es otro método, declarado explícitamente. Sus
primitivas exactas no tienen un error de cuadratura que deba convergerse. Deben
verificarse sus condiciones de rama, dominio, soportes y máximos de interpolantes
y conservarse los valores usados; esta nota no anticipa su resultado numérico.

La cota cubre únicamente **absorción saliente hacia `E > E_max`**. No limita
emisión entrante desde electrones externos arbitrarios, cortes inferiores,
resolución del interior del catálogo, error temporal ni una trayectoria distinta.
Tampoco reemplaza la admisión global de la etapa 2.
