# Modelo vigente: cierre neutral efectivo BCS/KWT — etapa 5 completada

**Admisión exploratoria del 12 de octubre de 2026.** Etapas 5.1, 5.2 y 5.3 completas. Certificado autoritativo: `/home/jdiaz/scratch/stage5_final_evidence_completed_20261012/admission.json`, SHA-256 `ed30c8525ecb3952e4042eca3b84e9ea0c2c78f41d7e0f17d880779b68747c4a`. [Estado y resultados](../../status.md); [alcance de etapa 5](stage5/README.md).

La API pysnspd.effective evoluciona conjuntamente condensado, cuasipartículas, fonones y circuito; obtiene el potencial por neutralidad cuasiestática del mismo estado. Es una hipótesis mesoscópica útil para exploración y calibración. El legado está congelado en v1.0.0. El intento microscópico de primer orden Moyal falló sus controles fotónicos y se conserva rechazado; no integra esta física efectiva. El [contrato anterior](../releases/history/MODELO_VIGENTE_before_v1.1.0_20261011.md) conserva la historia.

## Estado y energía

$$
X=\{d_i,p_{ik},n_{iq},I_b,I_s,v_c\},\quad d_i=\Delta_i/E_u,\quad E_u=k_BT_c,
\quad E_{ik}=E_u e_{ik},\quad e_{ik}=\sqrt{\bar\xi_k^2+|d_i|^2}.
$$

$i$ es nodo 2D; $k$ estado nativo electrónico fijo $\bar\xi_k>0$; $q$ modo fonónico. Las ocupaciones satisfacen $0\le p\le1$, $n\ge0$. $\phi_i$ [V] es algebraico. $V_i=d_fA_i^{dual}$ [m³], $u_i^*=N_0E_u^2V_i$ [J]; $N_0$ es DOS normal monoespín [J⁻¹m⁻³].

$$
U_{e+\nabla}=\sum_i u_i^*\left[\frac{|d_i|^2}{g}
-2\sum_kw_k(e_{ik}-\bar\xi_k)+4\sum_kw_ke_{ik}p_{ik}\right]
+\sum_{\langle ij\rangle}C_{ij}\kappa_{ij}|d_i-e^{-i\alpha_{ij}}d_j|^2,
$$
$$
g^{-1}=\sum_kw_k\frac{\tanh(\bar\xi_k/2)}{\bar\xi_k},\quad
C_{ij}=N_0E_u^2d_f\ell_0^2c_{ij},\quad
\kappa_{ij}=\tfrac12[\kappa(|d_i|^2)+\kappa(|d_j|^2)],
$$
$$
K(a;T_b)=\frac{\pi N_0\hbar D}{4a}\tanh\frac{a}{2k_BT_b},\qquad
\kappa(\rho)=K(E_u\sqrt\rho;T_b)/(N_0\ell_0^2),
$$
$$
{\cal F}_i=\partial_{d_i^*}U_{e+\nabla},\qquad
I^s_{ij}=-\frac{2e}{\hbar}\partial_{\alpha_{ij}}U_{e+\nabla}.
$$

$\omega_q=\Omega_q/E_u$ y $c_q=g_{ph}(\Omega_q)W_q^{ph}/(N_0E_u)$ son energía fonónica normalizada y capacidad nativa; $W_q^{ph}$ es el peso físico de cuadratura [J].

La energía BCS usa la misma cuadratura que fuerza y emparejamiento. La rigidez radial es **una extensión energética propia** motivada por corriente BCS del límite sucio: $T_b$ fijo, derivada $\kappa'$ completa y límite cartesiano finito en $\Delta=0$. No es un funcional Usadel exacto fuera del equilibrio. El depairing entra por el gradiente; $\Gamma_{\rm DOS}=0$. Los casos actuales usan $\alpha_{ij}=0$. La transferencia inicial alcanza $T_{E,eq}\simeq20$ K por encima de $T_c=8.65$ K: movilidad KWT y rigidez $\kappa(T_b)$ se adoptan como cierre fenomenológico calibrable, sin derivación microscópica válida para todo ese rango y sus distribuciones.

## Evolución conjunta

Definiendo $z_i=d_i/|d_i|$ ($z_i=1$ si $d_i=0$), $f_i=z_i^*{\cal F}_i/u_i^*$:
$$
D_td_i=\dot d_i+\frac{2ie}{\hbar}\phi_id_i
=-\frac{z_i}{A_i\tau_K}
\left[R_i^{-1}\operatorname{Re}f_i+iR_i\operatorname{Im}f_i\right],
$$
$$
\vartheta_i=\max(T_{E,i},T_b)/T_c,\quad
A_i=\sqrt{(1+\vartheta_i)/2},\quad\tau_K=\pi\hbar/(8k_BT_c),
$$
$$
\tau_{{\rm sc},i}^{-1}=\vartheta_i/\tau_{ee}(T_c)+\vartheta_i^3/\tau_{ep}(T_c),
\qquad R_i=\sqrt{1+4E_u^2\tau_{{\rm sc},i}^2|d_i|^2/\hbar^2},
$$
$$
\dot p={\cal C}_{e-ph}+{\cal T}_E+
\frac{f_{\rm FD}(E,T_E)-p}{\tau_{\rm BGK}}+{\cal H},\quad
\tau_{\rm BGK}=\tau_{ee}(T_c)T_c/\max(T_E,T_b),
$$
$$
\dot n={\cal C}_{ph-e}-[n-n_{\rm BE}(\Omega,T_b)]/\tau_{\rm esc}.
$$

Las colisiones compartidas de dispersión/recombinación conservan energía electrónica+fonónica y balance FD/BE. El transporte intercambia energía física igual entre vecinos, no difunde directamente las etiquetas $\xi$. BGK conserva energía electrónica. $T_E,T_{ph}$ son funcionales de **energía equivalente**, no estados térmicos impuestos durante la evolución. Los fonones son locales con escape; no hay difusión fonónica lateral.

El calor físico se calcula de la movilidad y del Joule:
$$
Q_{\Delta,i}=\frac{2u_i^*}{A_i\tau_K}
[R_i^{-1}(\operatorname{Re}f_i)^2+R_i(\operatorname{Im}f_i)^2]\ge0,\qquad
P_i=Q_{\Delta,i}+P_{J,i}.
$$
La fuente de ocupación $\mathcal H$ deposita exactamente esa potencia. No se añade calor igual a un residuo energético.

## Neutralidad y circuito

$B$ es incidencia orientada; $P$ prolonga potenciales reducidos a cortes equipotenciales de ancho completo, $B_g=P^{\mathsf T}B$:
$$
B_g[I^s+GB_g^{\mathsf T}\Phi]=s_gI_s,\quad
G_{ij}=c_{ij}/R_\square,\quad
V_{\rm dev}=\phi_L-\phi_R,\quad
P_J=\sum G_{ij}(\phi_i-\phi_j)^2\ge0.
$$
$s_g$ inyecta $+1$ a izquierda y $-1$ a derecha. Una referencia fija el gauge; cada cara reparte su Joule por mitades. La conductancia normal constante es hipótesis constitutiva; no se resuelve modo de carga $h_T$.

Los cortes son reservorios superconductores rígidos: amplitud DC fija, fases Josephson independientes, $D_t\Delta_{\rm cut}=0$, $Q_{\Delta,\rm cut}=0$ y poblaciones al baño. Su intercambio de calor $P_{\rm cuts}$ se registra explícitamente. Los bordes transversales son aislados. El control de longitud mide la influencia de esos cortes.

$$
L_b\dot I_b=V_b-R_bI_b-v_c-R_l(I_b-I_s),\qquad
L_{k,\rm ext}\dot I_s=v_c+R_l(I_b-I_s)-V_{\rm dev},
$$
$$
C_c\dot v_c=I_b-I_s,\quad V_{\rm out}=R_l(I_b-I_s),\qquad
L_{k,\rm ext}=96~{\rm nH}+(5~\mu{\rm m}-L)\frac{\hbar}{2e}
\left.\frac{dq}{dI}\right|_{\rm DC}.
$$
$L_{k,\rm ext}$ se calcula una vez desde el nuevo DC y permanece fijo. No hay $L_k$ lumped dinámico. $I_b$ o $I_s$ no positivos detienen el ensayo con checkpoint; no se recortan ni invierten corrientes.

## Referencias y parámetros libres

La referencia utiliza $T_c=8.65$ K, $T_b=0.9$ K, $D=5\times10^{-5}$ m²/s, $R_\square=608$ Ω, ancho 80 nm y **espesor físico $d_f=7$ nm**. La separación $h_{\rm mesh}$ es resolución numérica, no espesor. Reducirla no cambia por sí mismo la densidad de corriente. Al barrer espesor físico a material volumétrico fijo ($N_0,D$ constantes), debe conservarse $R_\square d_f$: el ejemplo 7→4,9 nm cambia 608→868,57 Ω. Mantener $R_\square$ al cambiar $d_f$ cambia el $N_0$ inferido y no compara el mismo material. $\ell_0=\sqrt{\hbar D/(2k_BT_c)}$ y $N_0=(2e^2DR_\square d_f)^{-1}$.

**$\tau_{ee}(T_c),\tau_{ep}(T_c)$ son calibradores libres**, con nomenclatura Allmaras (2020) y nominales 5 ps/24,7 ps. No son constantes medidas ni centro/límites del rango por calibrar. $\tau_{\rm esc}=15$ ps es otra hipótesis de interfaz.

El material acústico 3D efectivo usa $v_{\rm avg}=4912$ m/s, $\Theta_D=664$ K:
$$
g_{ph}(\Omega)=\frac{3\Omega^2}{2\pi^2\hbar^3v_{\rm avg}^3}
\mathbf1_{(0,k_B\Theta_D]},\quad
\alpha^2F(\Omega)=\frac{\hbar\Omega^2}{\pi\tau_0^{e-ph}(k_BT_c)^3},
\quad\tau_0^{e-ph}=\frac{720\zeta(5)}{\pi^2}\tau_{ep}(T_c).
$$
Variar $\tau_{ep}(T_c)$ modifica vinculadamente colisiones y movilidad KWT. nbn-a2f-ph.dat conserva procedencia/hash pero **no proporciona estas tasas**: su DOS absoluta no está identificada. No hubo cálculos Quantum ESPRESSO/DFPT propios.

## Fotón y reloj

Cada run prepara su propio DC compatible; no importa como prueba la admisión del backend Usadel histórico.
$$
E_{\rm ret}=\eta_{\rm ret}hc/\lambda,\quad
a_i=\mathbf1_{i\notin{\rm cuts}}e^{-|\mathbf r_i-\mathbf r_0|^2/(2\sigma^2)},\quad
E_{{\rm ret},i}=E_{\rm ret}\frac{V_ia_i}{\sum_jV_ja_j}.
$$

Se invierte energía BCS+fonónica actual para inicializar FD/BE comunes a $T_i^0$, manteniendo el gap DC. La deposición ocurre una sola vez. $\mathbf r_0$ está en el centro longitudinal; $\eta_{\rm ret}=0.667$ y $\sigma_{\rm ref}=7.549834435$ nm son referencias exploratorias ajustables.

$\sigma_{\rm ref}=\sqrt{2D(0.57~{\rm ps})}$ extrapola una escala de cascada; no es una medida del ancho para 775/1550 nm. **$t=0$ es la transferencia efectiva posterior a la cascada**. 0,57 ps no es una ventana térmica ni un retraso calculado.

## Balance y criterios de admisión

$$
E_{\rm tot}=U_{e+\nabla}+\sum_{iq}u_i^*c_q\omega_qn_{iq}
+\tfrac12(L_bI_b^2+L_{k,\rm ext}I_s^2+C_cv_c^2),
$$
$$
\dot E_{\rm tot}=V_bI_b-R_bI_b^2-R_l(I_b-I_s)^2-P_{\rm esc}-P_{\rm cuts}.
$$

El observador diferencia independientemente esta energía. El auditor exige:

- Defecto integrado ≤1 % de $E_{\rm ret}+\int|P_{\rm ext}|dt$.
- Residuo instantáneo ≤$10^{-20}$ W+$10^{-10}S_P$, donde $S_P=|I_sV_{\rm dev}|+|P_{\rm esc}|+|P_{\rm cuts}|+|V_bI_b-R_bI_b^2-R_l(I_b-I_s)^2|$.
- Control temporal: Vout y mínimo nodal/DC al 2 %; latencia dentro de $\max(2\%|\tau|,0.1~{\rm ps})$. Sólo este eje tiene ese piso absoluto, que nunca se suma al reloj físico. Malla y dominio: Vout, latencia, L2 del campo/DC y déficit/pico real al 5 %.
- Espectro exploratorio, sólo 775 nm y paquete 54→80 estados (32/6→64/8 electrónico, 96→144 fonones y 48→96 transporte): L2 del campo/DC y déficit central/pico real ≤5 %; mínimo nodal diagnóstico. Vout ≤5 % por decisión explícita **posterior a los datos** (SHA `bda09a6aaf095a175440dc77ac76faef8eaab98072fce2f5b5aa8984eb35efa6`); latencia conserva 2 % y energía 1 %. Los resultados anteriores 2 % y 3/2/2, incluido el FAIL eléctrico al 2 %, permanecen íntegros. La resolución adoptada es 54 estados exploratorios; 80 es control opcional de sensibilidad, sin prueba de convergencia al continuo.
- Soporte físico, corrientes globales positivas, fuentes/datos identificados y energía inicial integrada espacialmente.

Los cuatro refinamientos de convergencia corresponden al escenario 775 nm y a sus intervalos declarados. 1550 nm se presenta como trayectoria nominal exploratoria conservativa; su contraste temporal tardío cubre 29,62–31,62 ps, sin auditoría independiente de los cuatro ejes durante 0–200 ps.

El voltaje se compara con escala $\max(|V_{\rm out}|_{\rm pico},V_{\rm umbral})$ y se informa también el error respecto al pulso. Para paso se conserva el mínimo nodal del gap normalizado por su DC al 2 %. La política espectral exploratoria exige los dos observables robustos al 5 % y conserva el mínimo nodal como diagnóstico; los resultados originales al 2 % y 3/2/2 permanecen íntegros y visibles. Para malla y dominio, la política v2 compara la norma L2 del campo de magnitud en la ROI común120×80 nm (normalizada por gap DC) y el error de déficit central relativo al pico real de la curva de referencia, ambos al5 %. La coordenada x se centra en el depósito; y conserva su valor físico. En cada muestra guardada se interpola el campo complejo en espacio y se toma su magnitud; después se interpola ese observable escalar en tiempo. La comparación compleja temporal anterior se conserva y distingue como artefacto de análisis cuando depende de la fase global entre muestras; no se usa para inventar dinámica, winding ni eventos no guardados. El déficit usa el interpolante del déficit nodal positivo: es una aproximación espacial, no la integral exacta de la magnitud del interpolante complejo. El mínimo nodal y su criterio anterior permanecen visibles como diagnóstico, sin reclasificar resultados históricos. Sólo los eventos guardados confirmados dan latencia; interpolar curvas de análisis no crea eventos.

La tolerancia no lineal nominal es 2e-7, contrastada con 2e-9. La campaña histórica 775 nm conserva todos sus segmentos: nominal sin guarda local hasta 12,4 ps, 1e-5 hasta 22,435 ps y 2,5e-6 después; temporal sin guarda hasta 12,42 ps, 1,25e-6 hasta 22,2525 ps y 3,125e-7 a partir de esa ancla temporal ejecutada. Las anclas de controles deben situarse en [22,185;22,435] ps y el candidato nunca puede tener guarda más laxa que la referencia (temporal: razón 1/8); dentro de esa banda la transición numérica está acotada, sin fingir aislamiento exacto del eje. El protocolo exige declaración SHA, linaje íntegro e idénticas 13 fuentes físicas; no permite ignorar parámetros genéricamente. Las estadísticas de pasos se informan después de ambas anclas, sin presuponer razón dos en cada paso ni orden dos global. Los comandos fresh publicados fijan cap 0,02 ps, NL 2e-7 y guarda 2,5e-6 desde cero para ambos fotones; un futuro control temporal fresh requiere declaración propia con cap 0,01 ps y guarda 3,125e-7.

$$
|E_{k+1}-E_k-\Delta W_{k+1}|\le\varepsilon_E\left(E_{\gamma,\mathrm{ret}}+\sum_{j\le k}|\Delta W_j|+|\Delta W_{k+1}|\right).
$$

La guarda se reconstruye del historial guardado desde su ancla declarada; los tramos previos no se certifican retrospectivamente. El1 % de error global continúa siendo un criterio independiente.

Los antecedentes parciales rechazados permanecen históricos: mínimo nodal espectral 2,276239 % frente al 2 % y déficit/pico de malla 6,583534 % frente al 5 %. En la ventana final 0–36,6 ps, el contraste espectral conserva su FAIL de Vout 2,597424 % frente al 2 %, y pasa el margen del 5 % adoptado después de examinar los datos. El mínimo nodal de malla presenta 22,5717 %/DC y sigue diagnóstico; no se reclasifica como convergente. Los cuatro controles finales pasan sus criterios vigentes y la revisión manual SHA certifica formación y latencia con extremo crítico de 33 ps. Todas las colas posteriores a 36,6 ps están fuera de la convergencia de cuatro ejes; se conservan las auditorías individuales de originales completos, sin inventar una cota posterior.

El punto medio implícito acopla todo el sistema. La limitación continua de eventos compartidos conserva estequiometría y soporte; su dependencia del paso se verifica por refinamiento. No es clipping. La convergencia cinética usa norma L1 de energía sin cancelación más guarda de respuesta del gap; los residuos de gap/circuito y el balance independiente siguen activos.

## Observación y límites

La campaña admitida usa umbral **positivo de 0,1 mV**, confirmación a cruce+0,5 ps y final a **cruce+5 ps**. Sin cruce, el techo es 200 ps y la latencia queda sin asignar. La hotbelt es un diagnóstico geométrico separado: una componente conectada con $|\Delta_i|\le0.2|\Delta_{i,\rm DC}|$ toca ambos bordes transversales; no equivale al evento eléctrico. En 775 nm aparece por primera vez en un snapshot a 15,3 ps, con bracket 15,2–15,3 ps.

Esta versión no acredita jitter estadístico, retraso óptico absoluto, reset eléctrico total ni calibración universal de la película. La admisión se limita a las hipótesis, parámetros y ventanas del certificado final; 1550 nm no hereda los cuatro ejes de 775 nm.

Ecuaciones principales 2D: [Principal.pdf](../../Matematica/Principal.pdf); auxiliares evaluadas durante la evolución: [Auxiliar.pdf](../../Matematica/Auxiliar.pdf); precomputados: [Precomputo.pdf](../../Matematica/Precomputo.pdf); entradas, calibradores y dependencias críticas: [Parametros.pdf](../../Matematica/Parametros.pdf). Fuentes en `sandbox/stage5_math_documents/`. El [informe final](../../output/pdf/implementation/Informe_final_etapa5_modelo_efectivo_20261011.pdf) presenta los resultados físicos y el postproceso. Fuentes primarias desarrolladas allí: BCS (1957), KWT, [Allmaras (2020)](https://thesis.caltech.edu/13748/) y [Vodolazov (2017)](https://arxiv.org/pdf/1611.06060). La rigidez radial energética y la discretización conservativa se identifican como hipótesis/implementación propias.


## Escala de corriente y alcance temporal

La auditoría bulk efectiva da $I_{\mathrm{dep}}=30.367443\,\mu\mathrm{A}$ para7 nm. La corriente física preservada $I=21.297206\,\mu\mathrm{A}$ corresponde a $I/I_{\mathrm{dep}}=0.701317$ en este cierre, una fracción diferente de la del cierre Usadel histórico. Es una predicción interna del modelo efectivo, sin equivalencia experimental asignada.

La latencia parte del estado postcascada: $\tau_{\mathrm{rel}}=t_{\mathrm{cruce}}-t_0$. La cascada anterior y retardos del setup aportan offsets todavía sin calibrar. Ausencia de cruce se informa como censura de la hipótesis y ventana estudiadas; no valida ni refuta por sí sola el experimento. Winding CCW con $A=0$ identifica vórtices de los estados guardados, sin convertirlos en hotbelt ni detección eléctrica.

## Resultados y cobertura finales

El control nominal 775 nm cruza a 31,94131686165 ps y termina a 36,94131686165 ps; 1550 nm completa 200 ps sin cruce ni hotbelt transversal según el criterio guardado. Los defectos energéticos globales máximos son 0,663463 % y 0,000144050 %, respectivamente, normalizados por energía retenida más trabajo externo absoluto acumulado. La conservación no se presenta como ajuste experimental.

Los cuatro ejes se comparan en 0–36,6 ps. La malla alcanza L2/DC=4,817296 % (<5 %); el espectro alcanza L2/DC=3,085255 %, déficit/pico=1,971201 % y Vout=2,597424 % (márgenes 5 %), con diferencia de latencia 0,313528 ps (<2 % de la referencia). [Estado](../../status.md) contiene los cuatro resultados. Todas las confirmaciones ocurren antes de 33 ps; quedan al menos 3,640719 ps comparados después de la última confirmación. Las colas de los cinco originales 775 nm posteriores a 36,6 ps no tienen convergencia de cuatro ejes, aunque sus balances, soporte, eventos y signo de corriente se auditan hasta sus finales respectivos.

La sensibilidad de anchura 775 nm usa σ/√2 y √2σ en 0–20 ps, con igual energía retenida. La inicialización FD/BE a temperatura común cambia también el reparto electrónico/fonónico: no es sensibilidad geométrica pura. Menor o mayor temperatura inicial por sí sola no ordena la respuesta eléctrica. No se asigna latencia ni censura a 200 ps a esos controles de 20 ps.
