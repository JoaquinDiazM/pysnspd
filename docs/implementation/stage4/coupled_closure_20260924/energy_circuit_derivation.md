# Energía, trabajo y circuito del ensayo acoplado

24 de septiembre de 2026. Auditoría del contrato B–D y de las implementaciones
experimentales. Este documento separa las identidades ya implementables de la
identidad de trabajo que debe comprobar el nuevo cierre dinámico. No declara
aprobado un transiente que todavía no se ha ejecutado.

## 1. Qué está fijado y qué falta

El circuito completo de la memoria ya está cerrado e implementado en
`pysnspd/experimental/electrical_ports.py`. La deposición de calor B.41 también
es una elección constitutiva explícita: no corresponde sustituirla por una
temperatura ni volver a elegir su forma espectral. Los factores materiales de
las colisiones no se deducen de una prueba de conservación.

La tarea pendiente es conservar conjuntamente energía y carga cuando el
**espectro espacial** responde a un gap complejo móvil y las ocupaciones no
son uniformes. B.46 fue deducida para un catálogo local. Que la fuerza Keldysh
espacial no sea un gradiente a distribución local fija **no demuestra que no
exista energía interna**: la distribución cambia de coordenadas cuando se mueve
el espectro. Hay que incluir ese trabajo en la ecuación cinética y verificar
una energía calculada independientemente de la integral de potencias.

La fuente publicada más próxima revisada para la conexión de calibre es
[Vargunin y Silaev, PRB 96, 214507, apéndice A](https://arxiv.org/abs/1709.01444).
Conserva explícitamente enlaces temporales de calibre, el desplazamiento
electroquímico y fuentes proporcionales a la variación del gap. Su aplicación
de flujo de vórtices emplea distribuciones cuasiestacionarias; copiar sólo sus
ecuaciones finales no suministra almacenamiento temporal general. La revisión
de [Belzig et al.](https://arxiv.org/html/cond-mat/9812297v2) distingue las
ecuaciones espectrales, los dos modos cinéticos y la autoconsistencia.

## 2. Conteo espectral: una identidad útil del oráculo actual

Sean $d_i=\Delta_i/(k_BT_c)$, $m_i$ las áreas duales adimensionales y
$z=\eta-iE$, con $E$ en unidades de $k_BT_c$. La acción compleja ya usada por
`retarded_spatial_usadel.evaluate` es

$$
\mathcal A=\sum_i m_i\{2z(1-g_i)-d_i^*f_i-d_i\widetilde f_i\}
+\sum_{ij}c_{ij}\{(f_i-U_{ij}f_j)
(\widetilde f_i-U_{ij}^*\widetilde f_j)+(g_i-g_j)^2\}.
$$

Sobre su raíz estacionaria, la regla de la envolvente da

$$
\partial_z\mathcal A=2\sum_i m_i(1-g_i),\qquad
\partial_E\operatorname{Im}\mathcal A
=2\sum_i m_i(\operatorname{Re}g_i-1).
\tag{EC.1}
$$

La regla exige que la variación espectral de los reservorios sea nula o que
sus términos de reacción se incluyan. Por ejemplo, una condición de contorno
BCS que depende de $E$ aporta su derivada al diferenciar toda la acción; no
puede tratarse como fija y móvil a la vez.

EC.1 identifica el conteo **global** de estados en exceso respecto del metal
normal con $\operatorname{Im}\mathcal A/2$, más una constante de extremo.
También da identidades mixtas de DOS, fuerza y corriente mediante una nueva
derivada. Es una comprobación útil de escala, signo y contornos que utiliza el
mismo grafo; no requiere una simulación temporal larga.

El conteo global permite seguir ocupaciones de energía comunes a toda la
región. Para una distribución espacial arbitraria no basta reemplazar cada
catálogo local por $x_i(E)=\int^E\rho_i(E')dE'$. Ahora $x_i$ depende de todos
los $d_j$. Mantener cada $p_i(x_i)$ fija genera derivadas no locales, mientras
que el momento anómalo de la fuerza Keldysh pesa la distribución del nodo.
La conexión cinética que transporta estados y energía es parte del cierre,
no una identidad garantizada por el cambio de variable.

## 3. Energía independiente que debe ponerse a prueba

En la aproximación espectral adiabática y con simetría electrón–hueco, una
construcción candidata para el sector neutro es

$$
\overline U_e[d,h_L]=\overline U_{\rm eq}[d,T_b]
-2\sum_i m_i\int_0^\infty E\rho_i(E;d)
\,[h_{L,i}(E)-h_0(E)]\,dE,
\quad h_0=\tanh\frac{E}{2T_b/T_c}.
\tag{EC.2}
$$

La unidad de energía es $U_0=N_0(k_BT_c)^2d_{\rm film}\ell_0^2$,
con $N_0$ por espín y $\ell_0^2=\hbar D/(2k_BT_c)$.
El factor $-2$ procede de $h_L=1-2p$ y del factor $4N_0$ de B.46.
El estado térmico fuera de autoconsistencia se calcula con el **mismo**
funcional espectral: $U_{\rm eq}=F_{\rm eq}-T\partial_TF_{\rm eq}$ a campos
fijos, con corte y contraterminos consistentes. La energía libre $F$ sola no
es energía interna. Una derivada de una suma finita debe declarar si mantiene
fijo el corte físico o el número de frecuencias; sus errores no se cancelan
automáticamente al cambiar de cuadratura.

EC.2 recupera la energía de excitaciones de B.46 en el límite uniforme. Es un
**observable candidato**, todavía no una demostración del cierre espacial de
carga: la convención electroquímica utilizada para $h_T$, la neutralidad y
los términos de energía de orden cuadrático asociados deben derivarse en el
mismo calibre. No se añade un término arbitrario en $\phi h_T$ para hacer
cerrar el balance.

La derivada de EC.2 a $T_b$ fijo contiene

$$
\dot{\overline U}_e
=D_d\overline U_{\rm eq}[\dot d]
-2\sum_i m_i\int E\{\rho_i\dot h_{L,i}
+(h_{L,i}-h_0)D_d\rho_i[\dot d]\}\,dE.
\tag{EC.3}
$$

Es posible comprobar EC.3 por diferencias direccionales de estados, usando
las tangentes espectrales exactas o resolviendo de nuevo las raíces. La
comparación decisiva es su contracción con el **RHS cinético completo**,
incluido el trabajo del gap y del potencial. Definir una variable
`U += integral(power)` y comparar después esa variable con la misma integral
sería una identidad del registrador; no comprobaría EC.2 ni la física.

## 4. Calor y transferencias que se contabilizan una sola vez

Para una ley de condensado KWT con fuerza cartesiana coherente, la potencia
irreversible es $Q_\Delta=v^T\Gamma v\ge0$, donde $v$ representa $D_t\Delta$.
En el modelo energético original, el condensado pierde $Q_\Delta$ y las
ocupaciones reciben exactamente esa potencia. El trabajo reversible
$\mathbf j_s\cdot\mathbf E$ puede ser negativo; no se convierte en calor
mediante un valor absoluto.

Cuando la corriente disipativa se conserva como la ley óhmica heredada,
$P_J=\sigma_n|\mathbf E|^2$. Entonces B.41 distribuye
$P_{\rm heat}=P_J+Q_\Delta$:

$$
\mathcal H(x)=\frac{P_{\rm heat}\,E(x)f_*(E)[1-p(x)]}
{4N_0\int E^2f_*(E)[1-p(x)]\,dx},\qquad
4N_0\int E\mathcal H(x)dx=P_{\rm heat}.
\tag{EC.4}
$$

La normalización usa el soporte espectral y el volumen real de la celda. En
un código a energía fija se usa $dx=\rho\,dE$ y además la deriva de niveles;
no se trata el cambio de DOS como una creación de partículas. La misma
inyección debe aparecer en los valores de población almacenados, no sólo en
un acumulador de calor.

Si la cinética transversal pasa a calcular toda la corriente disipativa, no
se puede agregar además una corriente normal $\sigma_n\mathbf E$ que
represente los mismos portadores. Análogamente, si su ecuación energética ya
contiene calentamiento eléctrico microscópico, EC.4 no lo vuelve a inyectar.
Hay que identificar qué parte de la ley óhmica y de la movilidad KWT queda
como cierre residual del modelo mesoscópico. No se declara doble conteo sólo
porque exista una corrección espectral temporal: se debe señalar el término
disipativo concreto que coincide.

Las colisiones electrón–fonón intercambian igual energía con signo opuesto;
el escape fonónico es una pérdida al exterior. La normalización material NbN
pendiente impide atribuir tasas absolutas nuevas, pero admite un control con
colisiones desactivadas o con el escenario sintético ya registrado, cuya
naturaleza se mantenga visible.

## 5. Circuito de la memoria y puertos

El estado es $(I_b,I_s,v_c)$, con $I_s$ la corriente **total** de la rama
del detector. No representa sólo supercorriente. Las ecuaciones implementadas
son exactamente las de la adenda CM:

$$
\begin{aligned}
L_b\dot I_b&=V_b-R_bI_b-v_c-R_L(I_b-I_s),\\
L_{k,\rm ext}\dot I_s&=v_c+R_L(I_b-I_s)-V_{\rm dev},\\
C_c\dot v_c&=I_b-I_s,\\
V_{\rm out}&=R_L(I_b-I_s),\qquad V_d=v_c+V_{\rm out}.
\end{aligned}
\tag{EC.5}
$$

El puerto tiene orientación pasiva izquierda–derecha:
$V_{\rm dev}=\phi_L-\phi_R$, entrada de corriente $I_s$ por el extremo
izquierdo y salida por el derecho. Su potencia es $I_sV_{\rm dev}$, nunca
$I_bV_{\rm out}$. El condensador separa $V_d$ de $V_{\rm out}$.
Los valores heredados son $R_b=10$ kΩ, $L_b=1$ µH, $R_L=50$ Ω y
$C_c=100$ pF. El valor de $V_b$ del ensayo debe inicializar la misma rama:
$V_b=R_bI_{\rm ref}+V_{\rm dev,ref}$; 0,300 V no es compatible por
definición con cualquier corriente elegida.

La energía exterior y su identidad son

$$
U_c=\tfrac12L_bI_b^2+\tfrac12L_{k,\rm ext}I_s^2+\tfrac12C_cv_c^2,
\qquad
\dot U_c+I_sV_{\rm dev}+R_bI_b^2+R_L(I_b-I_s)^2=V_bI_b.
\tag{EC.6}
$$

El error integrado debe emplear valores de estados y cuadratura de potencia
independientes. Para no restar grandes potencias constantes de polarización,
se calcula también la diferencia contra una trayectoria del equilibrio
inicial con el mismo circuito; se declaran ambas definiciones.

La inductancia exterior se fija **antes** del transiente. En una rama
estacionaria estable, la diferencia de fase entre los mismos planos del puerto
define $\Phi_{\rm res}=(\hbar/2e)(\theta_R-\theta_L)$ y

$$
L_{k,\rm res}^{\rm diff}=\frac{d\Phi_{\rm res}}{dI},\qquad
L_{k,\rm ext}=10\ {\rm nH}-L_{k,\rm res}^{\rm diff}>0.
\tag{EC.7}
$$

La respuesta debe obtenerse de la rama espacial nueva; la partición de una
geometría o funcional anteriores no se trasplanta. No se actualiza EC.7 a cada
instante para forzar una inductancia total constante.

## 6. Contactos: un trabajo que no debe desaparecer

Con $D_t=\partial_t+2ie\phi/\hbar$, una fase de reservorio materialmente
estacionaria satisface $\dot\theta=-2e\phi/\hbar$. Por tanto

$$
\frac{d}{dt}(\theta_R-\theta_L)=\frac{2e}{\hbar}V_{\rm dev}.
\tag{EC.8}
$$

Fijar ambas fases en el laboratorio mientras los potenciales terminales son
distintos describe reservorios externamente actuados. Puede servir como
control, pero exige registrar su reacción y su potencia. Del mismo modo,
seguir una amplitud de reservorio $|\Delta_b(I_s)|$ introduce trabajo radial.
El trabajo de reservorios ya presente en $\mathbf J_{\rm tr}$ no desaparece
por utilizar una malla dual o por llamar al contacto «térmico».

El balance global que hay que recuperar es

$$
\frac{d}{dt}(U_e+U_{\rm ph}+U_c)
=V_bI_b-R_bI_b^2-R_L(I_b-I_s)^2-P_{\rm esc}
-\oint(\mathbf Q_e-\mathbf J_{\rm tr})\cdot\mathbf n\,dS.
\tag{EC.9}
$$

Aquí $Q_\Delta$ ya se canceló internamente y el trabajo del puerto se canceló
entre película y circuito. El signo de cada potencia debe salir de las
orientaciones declaradas, nunca elegirse posteriormente por valor absoluto.

## 7. Ensayo mínimo que mezcla la física

Se propone una cinta **2D** suave de 80 por 160 nm sobre la malla dual ya
admitida, sin fotón, con el material de referencia K20 y un equilibrio con
corriente no nula situado dentro de una rama estable. La longitud es la del
control anterior, no el rango final de confianza geométrica de Korzh.

El estado inicial combina una perturbación suave de amplitud, una distribución
electrónica con forma energética no térmica y una perturbación de fase con
variación transversal. Todos sus perfiles desaparecen en los contactos. Puede
iniciarse $h_T=0$, pero debe permitirse que la supercorriente y los gradientes
de $h_L$ lo generen. Al polarizar la rama, los términos cruzados aparecen ya
en su respuesta débil: no se superponen dos controles desacoplados a corriente
cero. La amplitud se registra antes de comparar pasos y mantiene el soporte
$|h_L\pm h_T|\le1$ sin recortes.

Una campaña económica tiene tres comparaciones realmente distintas:

1. **Estado de referencia:** conserva la rama, el punto fijo circuital y el
   balance de puerto; determina EC.7. No necesita esperar una recuperación
   de nanosegundos.
2. **Transiente mixto y mitad de paso:** guarda gap, ambas distribuciones,
   potencial, corriente del puerto, $V_{\rm out}$, calor que efectivamente
   entró a las poblaciones y EC.9. El paso KWT heredado sigue siendo de primer
   orden; el resto del acoplamiento no hereda automáticamente otro orden.
3. **Una variación de resolución dirigida:** refina la cuadratura energética
   o la malla que el piloto identifique como dominante. No repite todo el
   producto cartesiano de controles anteriores.

La eliminación instantánea de $h_T$ puede ejecutarse como comparación de
coste y efecto, identificada como reducción; la referencia conserva su
almacenamiento. Una comparación térmica/fotónica distinta no sustituye este
ensayo. Para las primeras trayectorias basta 1–5 ps si el piloto demuestra
intercambio medible; no se simula recuperación total para inflar el alcance.

Las observables aceptan inicialmente el margen práctico del 2 % de su escala
de señal, con suelo absoluto declarado cuando el observable tiende a cero.
El balance se normaliza por la energía perturbada o el trabajo realmente
intercambiado, no por la enorme energía inductiva de polarización. Una omisión
de un término de trabajo no se convierte en admisible relajando esa cifra.

Antes de consumir una hora, el control direccional EC.3 en una configuración
polarizada no uniforme debe revelar que la ecuación y la energía candidatas
se corresponden. Si no se corresponden, la salida útil es ese diagnóstico,
no otra trayectoria larga del cierre incompleto.

## Referencias internas y alcance

- B.41: forma espectral efectiva del calentamiento, ya seleccionada.
- B.43–B.51: conteo local adiabático, energía y fuerzas del modelo 0.4.
- C.14–C.21: movilidad KWT, calor y balance del funcional original.
- [Adenda del circuito](../../../modelo_v0_4/actualizaciones/circuito_memoria_20260922.md):
  topología CM, signos, planos de puerto e inductancia exterior.
- `longitudinal_reciprocal.py`: disponibilidad cuadrática del sector real sin
  corriente; no constituye un cálculo de $U_e$ ni del calor de EC.4.
- `frozen_kinetic_usadel.py`: corriente y torque compatibles con una identidad
  de carga estática; no suministra los términos temporales de EC.3.

Una identidad algebraica de circuito o de flujo no calibra tiempos materiales.
El pase del futuro ensayo acreditará el cierre y su implementación dentro de
la aproximación registrada; no acreditará la cascada fotónica ni la latencia
óptica de Korzh, que permanecen fuera de esta etapa sin fotón.

## Instrumentación incorporada y comprobaciones ligeras

`pysnspd/experimental/adiabatic_energy_observer.py` implementa EC.2–EC.3 en
unidades SI: recibe $U_{\rm eq}$ y su derivada del cálculo físico, y observa el
momento espectral sobre una cuadratura de energía fija. No fabrica esa energía
de referencia, no añade fuerza, no cambia poblaciones y no corrige el residuo
del balance. También observa energía fonónica, suma potencias externas de CM.9
y compara diferencias de energía guardadas contra trabajo por trapecios en
los pasos aceptados. Las contribuciones electroquímicas adicionales que exija
el cierre dinámico deben entrar desde su derivación física, no desde una
normalización automática del observador.

`pysnspd/experimental/thesis_circuit_harmonic.py` resuelve las tres ecuaciones
de la memoria para una impedancia $Z_{\rm dev}(\omega)$ **proporcionada por el
modelo de película**. Usa amplitudes máximas y la convención
$e^{+i\omega t}$. Devuelve corrientes, voltajes y potencia compleja, con
$\langle P\rangle=\operatorname{Re}(VI^*)/2$. La inductancia exterior puede
construirse mediante EC.7. Un puerto con parte real negativa no se recorta ni
se relabela como pasivo; su interpretación corresponde al punto de
polarización de la película. Esta utilidad no sustituye la película por una
resistencia escogida ni acredita su impedancia.

Cinco pruebas nuevas pasaron localmente en 0,115 s, con un hilo numérico:

- La energía de una población BCS no térmica seguida por su número de estado
  coincide con el observador a energía física fija. Su derivada coincide con
  el trabajo a conteo fijo y con diferencias independientes de estados.
- Omitir deliberadamente $\dot\rho$ produce un defecto detectado: no se
  compensa modificando $\dot h_L$ ni una fuente de calor.
- Rotar fases cambia los propagadores anómalos BCS pero mantiene la energía
  observada; omitir potencia externa en una historia también deja un residuo.
- La respuesta armónica reproduce el RHS temporal existente de CM.4 a tres
  frecuencias. Su balance de potencia se contrasta además con muestras
  temporales de un ciclo, incluyendo almacenamiento inductivo y capacitivo.
- Una partición inductiva no positiva se rechaza; una impedancia activa
  conserva su potencia negativa y cumple la misma identidad de circuito.

Las pruebas son controles de la instrumentación. No reemplazan la comprobación
del trabajo del oráculo espacial no térmico ni del cierre dinámico de carga.
