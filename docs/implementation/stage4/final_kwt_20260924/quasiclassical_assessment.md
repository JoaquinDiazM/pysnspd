# Qué puede hacerse instantáneo en la referencia NbN de Korzh

Investigación del 24 de septiembre de 2026. Alcance: decidir qué separación
de escalas respalda el ensamblaje dinámico de etapa 4; no calibrar otra
película ni validar todavía una detección fotónica.

**La descripción cuasiclásica difusiva es una reducción apropiada para
continuar este modelo efectivo. No convierte automáticamente la distribución
de energía ni el desequilibrio electrón–hueco en estados instantáneos.** La
ruta práctica recomendada conserva sus ecuaciones temporales y aprovecha el
espectro instantáneo sólo como aproximación adiabática explícita. Esa ruta
no exige empezar por un solver microscópico atómico ni por una convolución
completa de dos tiempos. Sí requiere terminar consistentemente los términos
de espectro móvil, fase y potencial antes de declarar cerrada la etapa 4.

## Cuatro reducciones diferentes

| Reducción | Qué se elimina | Criterio físico relevante | Qué no se elimina por ello |
|---|---|---|---|
| Cuasiclásica | Oscilaciones en la longitud de Fermi y dependencia radial rápida alrededor de la superficie de Fermi | Energías, frecuencias y gradientes de interés pequeños frente a las escalas de Fermi; metal suficientemente lejos de localización fuerte | Memoria superconductora y poblaciones fuera del equilibrio |
| Difusiva o Usadel | Anisotropía angular rápida del movimiento electrónico | Camino libre elástico mucho menor que la escala espacial; frecuencia por tiempo elástico mucho menor que uno | Difusión espacial, conversión de ramas y relajación inelástica |
| Espectro adiabático | Retardo del espectro respecto del condensado y los campos | Pequeñas correcciones temporales en las energías y lugares que contribuyen al observable | Trabajo reversible del espectro móvil ni evolución de sus ocupaciones |
| Distribución de carga instantánea | Memoria del desequilibrio electrón–hueco hT | Relajación de los modos de carga excitados mucho más rápida que su forzamiento | No se deduce de ninguna de las dos primeras filas |

[Belzig et al., §§2.2 y 2.4](https://arxiv.org/html/cond-mat/9812297v2)
presentan la reducción cuasiclásica y el límite sucio por separado; su
restricción expositiva a problemas estacionarios no restringe la teoría a
ellos. El texto de [Kopnin, capítulos 8–10](https://academic.oup.com/book/4832)
separa también teoría no estacionaria, aproximación cuasiclásica y ecuaciones
cinéticas. No se han importado fórmulas del libro que no estén accesibles en
las fuentes abiertas siguientes.

[Sauls, versión 2, §II, ecuaciones 7–12; §IV, ecuaciones 47–50](https://arxiv.org/pdf/2202.02260)
conserva el producto de convolución energía–tiempo tanto antes como después
de la reducción a Usadel. Su §IV.B elimina esa convolución para equilibrio;
la isotropización por impurezas no la elimina. Su §III, ecuaciones 27–31,
explicita los criterios difusivos. La expansión temporal de orden bajo es
una aproximación adicional. [Larkin y Ovchinnikov, página impresa 961,
alrededor de ecuaciones 8–9](https://www.jetp.ras.ru/cgi-bin/dn/e_041_05_0960.pdf)
introducen expresamente un régimen de baja frecuencia antes de simplificar
las funciones de dos tiempos.

## Datos utilizados y escalas calculadas

La referencia es la tabla suplementaria 1 de
[Korzh et al. 2020, PDF página 33](https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf):
ancho 80 nm, espesor 7 nm, Tc = 8,65 K, D = 0,5 cm²/s y resistencia de
hoja 608 Ω. La operación considerada es a 0,9 K. Son los parámetros de su
simulación; no todos son mediciones independientes. La misma tabla usa
τee(Tc) = 6 ps y τep(Tc) = 24,7 ps; el suplemento identifica ajustes del
modelo. **El jitter publicado de 2,6 ps es una dispersión de tiempos de
detección: no es la latencia media ni el tiempo local de variación del gap.**

Con la relación de gap débilmente acoplado **adoptada por nuestro modelo**,
Δ0 = 1,764 kB Tc, se obtiene:

| Cantidad calculada | Valor | Interpretación |
|---|---:|---|
| kB Tc | 0,74540 meV | Unidad energética de la referencia |
| Δ0 | 1,31488 meV | Gap de referencia, no nueva espectroscopia de la muestra |
| ℏ/Δ0 | 0,50059 ps | Escala coherente superconductora; no tiempo de colisión ni de carga |
| ℏ/(2 kB Tc) | 0,44152 ps | Unidad temporal usada en los controles espectrales |
| √(ℏ D/Δ0) | 5,00293 nm | Escala difusiva superconductora con esta convención |
| ℏ/(kB Tb), Tb = 0,9 K | 8,48693 ps | Escala asociada al ancho energético térmico; no tiempo de termalización |

La comparación algebraica para una frecuencia característica ω = 1/tvar
(tvar es tiempo de cambio, **no período sinusoidal**) es:

| tvar supuesto | ℏω/Δ0 | Lectura preliminar |
|---:|---:|---|
| 1 ps | 0,501 | No hay separación amplia para una variación de orden uno |
| 5 ps | 0,100 | Posible aproximación adiabática, dependiente del estado y observable |
| 10 ps | 0,0501 | Separación más favorable en la referencia con gap abierto |
| 50 ps | 0,0100 | Régimen claramente más lento respecto de esta escala |

Estos tiempos son escenarios para interpretar una trayectoria futura, no
latencias medidas ni nuevos límites impuestos al integrador.

## El límite difusivo está mejor justificado que la eliminación temporal

[Sidorova et al., tablas I–II y §IV](https://arxiv.org/pdf/1907.05039)
recopilan y miden películas NbN distintas: se encuentran tiempos elásticos
de 0,1 a 3,86 fs y caminos libres subnanométricos. Son comparadores de
material, no un intervalo de confianza para la muestra de Korzh. Con esos
tiempos y nuestro gap, Δ0 τel/ℏ va de 0,00020 a 0,00771; para variaciones de
1 ps, ωτel va de 0,00010 a 0,00386. Esto motiva eliminar la anisotropía
elástica incluso cuando se conservan procesos superconductores de
picosegundos. El ancho de 80 nm no es la única escala espacial: un núcleo
puede variar en aproximadamente 5 nm.

[Lomakin et al., tabla I y apéndice C](https://arxiv.org/pdf/2207.05012)
aportan otro conjunto, para películas de 2,5 nm: EF = 5–6,9 eV,
τel = 0,15–0,7 fs y kFℓ = 1,6–6,3, con parámetros electrónicos estimados
mediante transporte. Tomados por película, EF τel/ℏ es aproximadamente
1,6–6,3, no una jerarquía uniformemente muy grande. En cambio,
ℏ/EF = 0,095–0,132 fs es muy pequeño frente a picosegundos. Su velocidad
de Fermi ajustada de 5,3×10⁵ m/s, combinada **condicionalmente** con nuestro
D mediante τel = 3D/vF² y ℓ = 3D/vF, da 0,534 fs y 0,283 nm.

Estos comparadores apoyan separar las escalas atómicas de las
superconductoras. No sustituyen un conjunto de banda y tiempo elástico de
la muestra de Korzh ni acreditan EF τel/ℏ ≫ 1 en ella. Se conserva Usadel
como modelo efectivo del material, sin atribuirle exactitud microscópica
garantizada en NbN muy desordenado. No se usa la masa electrónica libre como
si fuese una masa de banda medida.

## Qué significa espectro suficientemente lento

Los [datos algebraicos reproducibles](scales/scales.json) acompañan la
[figura de jerarquía de escalas](scales/scale_hierarchy.png). La figura
compara escalas electrónicas de otras películas con la referencia Korzh y
grafica el indicador adiabático frente al tiempo relativo de cambio
|Δ|/|∂tΔ|, para diferentes amplitudes de gap. No grafica un transiente,
latencia ni jitter; tampoco estima τQ.

No basta comparar la duración total de una corrida con 0,50 ps. Para una
amplitud de gap no nula, un indicador local de cambio radial es

\[
 \epsilon_\Delta =
 \frac{\hbar\,|\partial_t|\Delta||}{|\Delta|^2}.
\]

Es una estimación dimensional de correcciones adiabáticas, no una cota
universal del error. Una variación del 1 % en 1 ps alrededor de Δ0 da
aproximadamente 0,005; una variación de orden uno en ese tiempo da
aproximadamente 0,5. **Un ensayo débil de 1 ps puede ser adiabático aunque
una supresión fuerte en 1 ps no lo sea.** El paso numérico elegido no altera
esta comparación física.

También importan la fase y el potencial en combinaciones invariantes de
calibre, el movimiento espacial y la anchura real de los rasgos espectrales.
Si un rasgo relevante cambia sobre una energía δE, la expansión temporal
puede contener ℏ/(δE tvar). El gap de referencia no controla uniformemente
un gap deprimido o una estructura espectral estrecha. Cerca de Δ = 0 el
indicador anterior pierde utilidad; no se lo regulariza para declarar
validez. Deben examinarse los coeficientes físicos que efectivamente entran
en los momentos.

La escala térmica de 8,49 ps advierte de derivadas energéticas estrechas,
pero **no impone esperar 8,49 ps a toda respuesta**: las energías de peso
despreciable en un espectro gapped no fijan por sí solas el error del
observable. Tampoco el ensanchamiento numérico η es un baño físico ni puede
usarse para inventar un tiempo de relajación.

[Vodolazov 2017, §II y ecuación 1](https://arxiv.org/pdf/1611.06060)
distingue la escala ℏ/|Δ| de la cascada y conserva una ecuación temporal de
ocupación con trabajo por cambio de amplitud. Sus simplificaciones para la
cascada y los modelos espaciales posteriores no autorizan omitir ese
trabajo en una distribución espacial no térmica general.

## La relajación de carga no es el tiempo elástico

hL describe el modo de energía; hT, el desequilibrio entre ramas de tipo
electrón y hueco. Su conversión y transporte dependen de energía, gap,
corriente, contactos y colisiones. Un tiempo elástico corto justifica
difusión, no la condición ∂t hT = 0. La neutralidad eléctrica tampoco
equivale a hT = 0: el potencial y el momento de distribución intervienen
conjuntamente, como muestran
[Golub, ecuaciones 9–10](https://www.jetp.ras.ru/cgi-bin/dn/e_044_01_0178.pdf).
Su límite próximo a Tc no se usa aquí como calibración a 0,9 K.

[Hübler et al., §II, ecuaciones 3 y 8–9](https://arxiv.org/pdf/1002.0983)
describen un tiempo de desequilibrio dependiente de energía y mecanismos
distintos de conversión. Sus medidas son en aluminio y no se trasladan a
NbN. Para NbN, [Sofer et al., discusión del modelo de magnetorresistencia](https://d-nb.info/1282362720/34)
también relacionan conversión, gap, campo y relajación inelástica; sus
nanohilos y régimen experimental tampoco fijan τQ para Korzh.

No se ha encontrado τQ(E, estado) de esa muestra que demuestre
τQ/tvar ≪ 1 en toda la ventana de pocos picosegundos. Las tasas τee y τep
de la tabla de Korzh no son τQ. Sus valores a Tc tampoco son ya mucho
menores que 1 ps. No se extrapolan a 0,9 K como mediciones ni se sustituye
τQ por el tiempo de apantallamiento eléctrico.

Como ilustración geométrica propia, la difusión normal de un modo seno
entre contactos absorbentes separados 80 nm tiene L²/(π²D) = 12,97 ps.
**No es una estimación del τQ superconductivo:** la conversión puede
acelerarlo y los coeficientes espectrales pueden cambiarlo. Sólo demuestra
que colisiones elásticas de femtosegundos no hacen instantánea una
población espacial.

## Recomendación para seguir sin añadir una cadena de pruebas

1. **Conservar las reducciones ya útiles:** movimiento cuasiclásico,
   isotropización difusiva y malla dual. No resolver oscilaciones de Fermi ni
   la cascada óptica excluida del alcance. Las energías retenidas deben
   permanecer en la banda de validez de esa reducción.
2. **Conservar temporalmente hL y hT.** Esta elección evita inventar una
   separación τQ ≪ tvar. Si una parte es numéricamente rígida, tratarla de
   forma implícita mantiene el sistema físico; eliminarla algebraicamente
   exige la separación adicional. No es necesario convertir un ensayo de
   tasa desconocida en otra condición previa al desarrollo.
3. **Reutilizar el oráculo espectral como cierre adiabático declarado** y
   derivar la cinética de primer orden temporal, incluyendo trabajo,
   neutralidad y fase/potencial. Una extensión completa de dos tiempos se
   reserva para un régimen en que esta aproximación falle de forma material
   para la señal, no se impone ahora. Suponer que la aproximación cubrirá
   toda la formación del hotbelt sigue pendiente de una trayectoria física.
4. **El próximo ensayo debe probar esa unión**, con campos suaves y una
   polarización definida; no otra suma de controles desacoplados. Mantener
   el circuito completo, el avance KWT admitido y la política práctica de
   error. Relajar una tolerancia numérica no sustituye una identidad física
   de trabajo o carga.

La recomendación distingue evidencia de hipótesis: la isotropización
rápida tiene respaldo material; la ausencia de memoria de hT a pocos
picosegundos no lo tiene. El espectro adiabático es una aproximación útil y
condicional, cuya extensión general aún exige el ensamblaje descrito en
[la revisión del cierre](coupled_closure_audit.md). Esta investigación no
resuelve por sí sola el trabajo espectral de B.46–B.51 trasladado del
catálogo local al oráculo espacial, ni los términos de calibre. No se ha
modificado la física ejecutable ni declarado cerrada la etapa 4.
