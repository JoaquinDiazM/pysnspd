# Núcleo real: correspondencia con el límite Ginzburg–Landau

Se calculó un contraste nuevo del perfil longitudinal real que atraviesa
`Delta=0`, a corriente nula y temperatura impuesta. La referencia es la
solución analítica `tanh` del funcional Ginzburg–Landau (GL), relacionada con
el límite de corriente nula de Langer–Ambegaokar. El candidato usa el potencial
térmico BCS completo y el gradiente local `K0`. El cálculo es independiente
de la malla GLL, del catálogo retardado y de la cuadratura cinética del solver.

**Resultado:** al acercarse a `Tc`, la energía excedente y la fuerza del
perfil candidato se acercan al límite GL conocido. Esta es una comprobación
de correspondencia física en un régimen controlado; no valida el núcleo de
Korzh a 0.9 K.

| Temperatura `T/Tc` | Exceso de amplitud BCS respecto de GL | Diferencia de energía BCS/GL | Máxima fuerza residual escalada muestreada del perfil BCS |
|---:|---:|---:|---:|
| 0.95 | 1.783 % | 2.750 % | 0.652 % |
| 0.98 | 0.694 % | 1.066 % | 0.256 % |
| 0.99 | 0.344 % | 0.527 % | 0.127 % |

Las temperaturas son controles independientes de correspondencia. El caso
de dispositivo conserva `T=0.9 K`, `Tc=8.65 K` y todos sus parámetros. Tampoco
se adopta una reducción espacial 1D para la cinta de 80 nm: la referencia
longitudinal se usa solamente como un problema cuyo perfil se conoce.

## Referencia y candidato

Se define `t=T/Tc`, la coordenada `x` en unidades de `ell0` y el campo real
con signo `d=Delta/(kB*Tc)`. La energía longitudinal se mide en
`N0*(kB*Tc)^2*A*ell0`, donde `A` es el área transversal. Para el potencial
térmico BCS homogéneo se evalúa

\[
V_{\rm BCS}(d)=d^2\ln t+2\pi t\sum_{n\ge0}
\left[\frac{d^2}{\epsilon_n}-2\left(\sqrt{\epsilon_n^2+d^2}-\epsilon_n\right)\right],
\qquad \epsilon_n=2\pi t(n+1/2).
\]

El término de la suma se calcula como
`d^4/[epsilon_n*(sqrt(epsilon_n^2+d^2)+epsilon_n)^2]`, expresión idéntica que
evita restar números casi iguales. El gap lejano `d_BCS` minimiza este
potencial y se obtiene resolviendo su ecuación autoconsistente.

La expansión cuártica define el control GL:

\[
V_{\rm GL}(d)=-\alpha d^2+\beta d^4,
\qquad \alpha=-\ln t,
\qquad \beta=\frac{7\zeta(3)}{16\pi^2t^2},
\qquad \kappa=\frac{\pi}{4}.
\]

Se retienen `ln(t)` y el factor `1/t²` de la expansión térmica; en el límite
crítico se recuperan `alpha≈1-t` y el coeficiente cuártico evaluado en `Tc`.
Ambas energías usan el mismo gradiente `kappa*(d')²`. Este acuerdo no es una
validación independiente de `K0` fuera del régimen GL.

La ecuación estacionaria GL es
`-2*kappa*d''-2*alpha*d+4*beta*d³=0`. Por sustitución directa,

\[
d_{\rm GL}(x)=d_{{\rm GL},\infty}
\tanh\frac{x}{\sqrt2\,\xi_{\rm GL}},
\qquad d_{{\rm GL},\infty}=\sqrt{\frac{\alpha}{2\beta}},
\qquad \xi_{\rm GL}=\sqrt{\frac{\kappa}{\alpha}}
\]

la satisface en toda la recta. Su energía excedente exacta es

\[
\Delta F_{\rm GL}=\frac{8\sqrt2}{3}
\frac{\alpha^2}{4\beta}\,\xi_{\rm GL}.
\]

La expresión se comprueba integrando `sech⁴`, cuya integral sobre la recta
es `4/3`. Esta derivación fija la normalización sin copiar una energía con
un convenio distinto de densidad de estados.

Para el candidato se prescribe
`d_BCS(x)=d_BCS,infinito*tanh[x/(sqrt(2)*xi_GL)]`: se conserva la anchura de
referencia y cada teoría usa su propio estado uniforme autoconsistente en
los extremos. Se resta el potencial de ese mismo estado para calcular la
energía excedente. Por tanto, la diferencia de energía incluye las
correcciones térmicas de orden superior y la pequeña diferencia de amplitud
lejana; no se interpreta como comparación de dos barreras estacionarias
exactas a condiciones lejanas arbitrariamente idénticas.

El perfil BCS está prescrito, no resuelto como punto de silla. Su fuerza
residual es `V'_BCS(d)-2*kappa*d''`; se divide por `2*alpha*d_BCS,infinito`
para obtener la última columna. La fuerza GL de su perfil analítico es
cero. El JSON contiene ambos perfiles y la fuerza para graficarlos.

## Qué aprende esta prueba y qué queda abierto

El límite térmico de un núcleo real queda conectado con una solución
analítica conocida, sin depender de que la malla experimental reproduzca
sus propias derivadas. Como el campo es real, `Im(conj(d)*d')=0` y
`q_delta=0`, incluso al atravesar el cero. La prueba no examina la
regularización de una fase que gira alrededor de un vórtice, ni decide el
valor físico de `delta`.

No hay evolución temporal, movilidad KWT calibrada, reparto de calor,
transporte no térmico, barrera del detector ni evento de fase observado.
El siguiente contraste físico pertinente debe estudiar la respuesta
espacial espectral en el régimen de baja temperatura y gradiente de fase
que esta prueba excluye. Las identidades y el refinamiento numérico se
mantienen como controles diferentes de esa comparación física.

## Coste y comprobación

El cálculo tomó 0.266 s, además de la carga de bibliotecas. Una repetición
con 256/512 frecuencias de Matsubara y 160/320 puntos de cuadratura cambió
la energía BCS menos de `3.9e-14` relativo. Las sumas incluyen sus colas
analíticas mediante zeta de Hurwitz. La cuadratura GL sobre ocho anchuras
a cada lado difiere de la integral de toda la recta en `1.2e-13`. Son
comprobaciones numéricas del control, no tolerancias exigidas a futuros
transientes ni incertidumbres materiales.

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/gl_core_reference.py \
  --output tmp/stage4_gl_core_reference_recheck.json
```

Fuentes primarias para el contexto y el funcional GL:
[Langer y Ambegaokar (1967)](https://doi.org/10.1103/PhysRev.164.498) y
[Skvortsov y Polkin, ecuación 4 y discusión del límite GL](https://arxiv.org/html/2506.18130v1).
La sustitución `tanh`, su normalización y la comparación BCS anterior se
desarrollan explícitamente aquí; no se atribuye a esas fuentes una validación
del candidato de baja temperatura.
