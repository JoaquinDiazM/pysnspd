# Puente longitudinal débil: capacidad preparada, ensayo temporal pendiente

El prototipo acopla una variación espacial real del gap con una distribución
electrónica que puede tener forma no térmica. El primer ensayo se restringe
a una referencia de fase constante y corriente nula: allí los bloques de
mezcla entre carga y energía desaparecen por simetría. No se impone esa
propiedad al vórtice sobre malla cartesiana ni al detector polarizado.

Se han ejecutado seis controles pequeños; todavía no se ha integrado esta
dinámica acoplada. La identidad que sigue verifica un bloque del modelo
mesoscópico débil, no cierra la etapa 4 ni valida el dispositivo.

## Variables, ecuaciones y origen del término cruzado

Se mantienen las unidades del núcleo térmico: energía en kBTc, área en
ℓ0² y tiempo τ=t/tD, con tD=ℏ/(2kBTc). Sean x=δ|Δ|/(kBTc), H el Hessiano
térmico de todos los nodos y μ el mapa inverso de movilidad radial KWT,
aplicado a fuerzas integradas en las áreas m_i. La referencia física del
primer ensayo debe satisfacer la ecuación del gap del mismo sumatorio
finito; el prototipo no borra el G0 no nulo del vórtice anterior.

Para cada energía de cuadratura E_j, con peso positivo w_j, se define

\[
h_0=\tanh\frac{E}{2t_b},\qquad
\chi=\partial_Eh_0,\qquad y=\delta h_L/\chi,\qquad
M_j=\operatorname{diag}(m_i\rho_{ij}).
\]

t_b=Tb/Tc y ρ=Re g son adimensionales. y tiene dimensión de energía en
estas unidades; es una coordenada de perturbación longitudinal, no una
temperatura ni un potencial electrostático. Su dependencia energética es
libre. Esta elección evita dividir numéricamente por una susceptibilidad
exponencialmente pequeña al formar la dinámica.

El kernel k_j se extrae de la misma fuerza Keldysh ya implementada:
δG=Σ_j w_j k_j δh_L,j. En el sector de fase constante, f=\tilde f y

\[
k_{ij}=\operatorname{Re}\{i m_i(f_{ij}-\tilde f_{ij}^{*})\}
=-2m_i\operatorname{Im}f_{ij}.
\]

Con L_j el bloque longitudinal del operador de caras existente, las
ecuaciones del prototipo son

\[
\dot x=-\mu\left(Hx+\sum_jw_j\chi_jk_jy_j\right),\qquad
M_j\dot y_j=L_jy_j+\tfrac12k_j\dot x.
\tag{L1}
\]

El producto k_j·dot x es nodo a nodo. En la coordenada δh_L, la segunda
ecuación se escribe M_j δdot h_L,j−(χ_j/2)k_j dot x=L_j δh_L,j.
Sustituir k=−2mR2 recupera el término ρ δdot h_L+R2 χ dot x que aparece
en la estructura cinética de [Vodolazov, ecuaciones 1 y 9](https://arxiv.org/pdf/1611.06060).
La ecuación 9 incluye gradientes espectrales. Esa correspondencia fija el
signo del término de amplitud; no demuestra que L1 incluya todos los
términos temporales de la teoría Keldysh general.

El factor 1/2 no es un ajuste: se obtiene al expresar el término R2 de la
ecuación 1 de Vodolazov mediante el kernel ya normalizado k=−2mR2 y usar
hL=1−2n. Adoptar L1 fuera de este sector débil es una elección constitutiva
que no se hace aquí. La cancelación de L3 y los controles discretos sí son
demostraciones algebraicas para el prototipo definido.

El coeficiente de L1 también se obtiene por reciprocidad energética.
La disponibilidad cuadrática, es decir, la expansión de la energía libre
cerca de la referencia térmica, es

\[
\mathcal A_2=\tfrac12 x^THx+
\sum_jw_j\chi_j y_j^TM_jy_j.
\tag{L2}
\]

La contribución de las distribuciones se obtiene de la segunda variación
de la entropía de Fermi: δ²(U−TbS)/2=mρ(δh_L)²/χ. Al derivar L2 y usar
L1, el trabajo cruzado se convierte exactamente en la fuerza no térmica
que ya recibe el gap:

\[
\dot{\mathcal A}_2=
\left(Hx+\sum_jw_j\chi_jk_jy_j\right)^T\dot x
+2\sum_jw_j\chi_jy_j^TL_jy_j
=-Q_{\rm KWT}-Q_L.
\tag{L3}
\]

En una referencia admitida, L_j posee forma de Dirichlet,
y^TL_jy=−Σ_aristas c_L,j(y_b−y_a)², con c_L,j≥0. Los incrementos de
gap y distribución se anulan en los contactos. El prototipo comprueba esa
estructura, los bloques de carga nulos y la positividad de la DOS; no
recorta coeficientes ni introduce capacidad donde ρ=0. Si la masa espectral
es cero, entrega la ecuación en forma de masa y rechaza dividir por ella.
Los operadores deben compartir el mismo gap real y declarar z=η−iE para
cada nodo energético. Se verifica su estacionariedad con esa energía,
además de exigir respuestas reales, finitas y de forma correcta al Hessiano
y a la movilidad. No se descarta silenciosamente una componente compleja.

## Qué conserva y qué falta

L3 es un balance de **disponibilidad**, no de energía interna total ni de
latencia. La positividad de las pérdidas no demuestra por sí sola estabilidad
de cualquier referencia: también importa el Hessiano H. Q_KWT es de orden
cuadrático en la perturbación. Devolverlo a las poblaciones requiere su
fuente de calentamiento de ese mismo orden o un balance explícito del baño;
no se añade dos veces ni se lo inventa como término de primer orden.

El presente módulo omite colisiones y fonones para aislar el intercambio
reversible gap–población y el transporte longitudinal. Ese límite se declara
en el plan, no sustituye las tasas ni los estados del futuro dispositivo.
La η usada para evaluar espectros es un regulador causal; no entra como
colisión, escape ni tiempo de relajación.

La coordenada local de conteo x_c(E,d)=∫₀ᴱρ(E',d)dE' puede conservarse en
una implementación posterior, pero requiere transformar el RHS. Si la
cinética física es ρ ∂tf|E+W ∂Ef=C+transporte y
v_c=−(∂t x_c)/ρ, entonces

\[
\rho\,\partial_tp|_{x_c}=C+\text{transporte}
+(\rho v_c-W)\partial_Ef.
\]

El último término sólo desaparece cuando la dinámica física lo justifica;
una DOS espacial no permite identificar automáticamente sus cuantiles con
niveles locales materiales. La formulación matricial del trabajo del gap
puede contrastarse con [Kozorezov y colaboradores, ecuaciones 9–13](https://link.aps.org/accepted/10.1103/PhysRevB.92.064504).
No se extrapola su reducción térmica/local a un remapeo espacial arbitrario.

## Archivos y siguiente ensayo

- Módulo: `pysnspd/experimental/longitudinal_reciprocal.py`.
- Controles: `tests/test_longitudinal_reciprocal.py`.
- [Plan mínimo](plan.json), [registro de pruebas](unit_tests.log) y
  [recibo con hashes](unit_tests_receipt.json).

El siguiente ensayo debe usar una referencia uniforme real, una perturbación
espacial de amplitud y otra de distribución con forma energética no térmica,
conservar todos los nodos, y observar el intercambio entre ambos campos.
El operador térmico completo aporta H; el kernel Keldysh aporta k y L.
Se comprobarán respuesta temporal, almacenamiento y balance L3 antes de
añadir colisiones, calor de segundo orden o corriente aplicada. No se
solicita otra ronda de precisión puntual de la DOS ni se presenta esta
preparación como una ejecución pendiente del usuario: aún falta el driver.
