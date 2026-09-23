# Referencia espacial lineal: sesgo del cierre local de gradiente

Se ejecutó una referencia nueva e independiente del catálogo y de los
operadores experimentales. Utiliza las ecuaciones C.29–C.31, con la solución
térmica BCS uniforme para `Tc=8.65 K`, `T=0.9 K`, `D=0.5 cm²/s` y corriente de
fondo nula. El cálculo tomó 0.032 s, más la carga del intérprete y bibliotecas;
no hubo pasos temporales ni transientes físicos.

**Resultado:** el cierre local de gradiente endurece las modulaciones cortas
respecto de la referencia Usadel espacial linealizada. El sesgo persiste en
el continuo; refinar la malla no lo elimina.

| Número de onda `k*ell0` | Longitud de onda (nm) | Rigidez Usadel `H_U` | Rigidez local `H_K0` | Exceso local | Tiempo local / Usadel, movilidad idéntica |
|---:|---:|---:|---:|---:|---:|
| 0 | Uniforme | 1.999999 | 1.999999 | 0 % | 1 |
| 0.15 | 196.81 | 2.019929 | 2.035342 | 0.763 % | 0.9924 |
| 0.30 | 98.40 | 2.078462 | 2.141371 | 3.027 % | 0.9706 |
| 0.50 | 59.04 | 2.210262 | 2.392698 | 8.254 % | 0.9238 |
| 1.00 | 29.52 | 2.728286 | 3.570795 | 30.881 % | 0.7641 |

`H` mide el cambio de la fuerza normalizada por unidad de modulación de la
amplitud. La última columna es un cociente de tiempos de relajación lineal,
no una latencia del detector. Se obtiene únicamente si se aplica la misma
movilidad positiva a las dos respuestas; no calibra esa movilidad. Un
cociente 0.764 indica una relajación aproximadamente 24 % más rápida bajo
ese supuesto, no un error demostrado de 24 % en el tiempo de detección.

La escala usada aquí es `ell0=sqrt(hbar*D/(2*kB*Tc))=4.69849 nm`. Las longitudes
de onda de la tabla son `2*pi/k`, no `1/k`. No se reutiliza la escala numérica
de las figuras históricas de C. Antes de interpretar longitudes cortas como
un margen físico del material debe comprobarse también `k*l << 1`, donde `l`
es el libre recorrido elástico; esta entrega no dispone de una medida que
cierre ese margen.

## Restricciones del estado y normalización

La temperatura está impuesta y el espectro se relaja térmicamente en la
referencia. No se mantienen congeladas las ocupaciones electrónicas. Por
ello no se equipara este resultado con la Hessiana a poblaciones fijas de
los ensayos 4A, ni se simula la cinética de B.

Con `t=T/Tc`, `d=|Delta|/(kB*Tc)`,
`epsilon_n=2*pi*t*(n+1/2)` y `E_n=sqrt(epsilon_n²+d²)`, el gap se obtiene de

\[
\ln t+2\pi t\sum_{n\ge0}
\left(\frac{1}{\epsilon_n}-\frac{1}{E_n}\right)=0.
\]

Se usa la convención BCS débil, `Delta0/(kB*Tc)=pi*exp(-EulerGamma)`. La fuerza
térmica uniforme es dos veces `d` por el miembro izquierdo; en su raíz su
derivada da

\[
H_U(0)=4\pi t\sum_{n\ge0}\frac{d^2}{E_n^3}.
\]

El incremento espacial es la primera expresión de C.30; el cierre local
añade `pi*(k*ell0)²/2` al mismo `H_U(0)`. Esto evita comparar estados
homogéneos o normalizaciones diferentes. La perturbación es real, pequeña
y sin gradiente de fase: `q_delta=0`. Por tanto el experimento no decide la
regularización del núcleo ni la nucleación de un vórtice.

## Comprobación numérica y reproducción

Se sumaron 2500 y 5000 frecuencias de Matsubara. La cola de alta frecuencia
se añadió con potencias inversas y la zeta de Hurwitz: órdenes 3, 5 y 7 para
el gap y `H_U(0)`, y órdenes 2 a 6 para el incremento espacial. Esta última
cola empieza en `1/n²`, y omitirla deja un error acumulado que decrece como
el inverso del corte.

El JSON conserva las sumas truncadas y las correcciones por separado. La
mayor diferencia relativa entre rigideces corregidas a ambos cortes es
`2.20e-16`. Ese número describe estabilidad numérica dentro de esta fórmula,
no una incertidumbre física ni un criterio exigido para los transientes.
La mejora obtenida con la corrección de cola evita convertir una suma más
grande en un cálculo pesado.

Script: `sandbox/stage4_core/linear_usadel_reference.py`.
Datos: `linear_usadel_reference.json` en esta carpeta. Para repetir sin
sobrescribir evidencia:

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  sandbox/stage4_core/linear_usadel_reference.py \
  --output tmp/stage4_linear_usadel_reference_recheck.json
```

Este diagnóstico proporciona una comparación física para la etapa 4. No
cierra por sí solo D.4.4: quedan la física del núcleo, el régimen no térmico,
la movilidad y el reparto de calor. Su uso inmediato es separar un error
de representación espacial de una diferencia constitutiva del modelo.
