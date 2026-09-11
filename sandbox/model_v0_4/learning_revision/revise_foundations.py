"""Apply the authored E01/E03 revisions; preserve the prior notebook first."""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[3]
DOCS = ROOT / 'docs/modelo_v0_4'
WORK = Path(__file__).resolve().parent
SOURCE = DOCS / 'E_cuaderno_de_aprendizaje_v0_4.md'
old = SOURCE.read_text(encoding='utf-8')
archive = DOCS / 'historial/E-r01'
archive.mkdir(parents=True, exist_ok=True)
for src in [SOURCE, DOCS/'latex'/f'{SOURCE.stem}.tex',
            ROOT/'output/pdf/modelo_v0_4'/f'{SOURCE.stem}.pdf',
            DOCS/'verificaciones/QA_visual_E.json']:
    dst = archive/src.name
    if src.exists() and not dst.exists():
        shutil.copy2(src, dst)

def section(start, end):
    return old[old.index(start):old.index(end)].strip()

e01 = section('# E.1.', '# E.2.')
start = e01.index('En mecánica,')
end = e01.index('El Modelo Estándar')
e01 = e01[:start] + r'''Un **sistema** es la parte del mundo que decidimos describir: una masa sujeta a un resorte, una cavidad electromagnética o un tramo de conductor conectado a una fuente. Elegimos sus límites, las interacciones relevantes y las variables que responden a nuestra pregunta. En mecánica, una coordenada $q(t)$ puede indicar la **posición** de una masa; su desplazamiento respecto de una posición de referencia $q_0$ es $u(t)=q(t)-q_0$. Otras coordenadas posibles son un ángulo o la elongación de un resorte. Coordenada y desplazamiento sólo coinciden si se elige así la referencia.

Un **campo** asigna un valor a cada punto de un dominio. Por ejemplo, $T(\mathbf r,t)$ asigna una temperatura a cada posición $\mathbf r=(x,y,z)$ y a cada instante $t$. El valor puede ser un número, como $T$, o varias componentes: el desplazamiento elástico $\mathbf u(\mathbf r,t)=(u_x,u_y,u_z)$ asigna un vector a cada punto. Las componentes son los valores que toma el campo; $\mathbf r$ y $t$ son las coordenadas donde se lo evalúa. El campo eléctrico $\mathbf E(\mathbf r,t)$ también es vectorial y no representa por eso la trayectoria de una partícula.

Una **configuración** es una elección de todos esos valores a un instante: el perfil completo de temperatura, por ejemplo. En una descripción clásica, el **estado** reúne los datos que el modelo necesita para continuar la evolución. Puede ser la configuración de un campo o de varios campos definidos en un mismo dominio, si éstos bastan. Para una cuerda vibrante se necesitan tanto $u(x,t_0)$ como su velocidad $\partial_tu(x,t_0)$; dos cuerdas con el mismo perfil y velocidades opuestas evolucionarán de manera distinta. Tampoco es obligatorio que todas las variables compartan dominio: un dispositivo puede combinar campos en un conductor con corrientes y tensiones de circuito que sólo dependen del tiempo. En una descripción cuántica, el estado incluye además la información probabilística y las correlaciones: no se reduce en general a asignar números clásicos a cada punto.

Un **modo** es una forma independiente de variación admitida por el modelo y sus condiciones de borde. En una cuerda fija en $x=0,L$, una forma posible es $\phi_1(x)=\sin(\pi x/L)$; su amplitud temporal $a_1(t)$ determina cuánto participa en $u(x,t)=a_1(t)\phi_1(x)$. Una **excitación** es una desviación respecto de un estado de referencia. Clásicamente puede consistir en aumentar esa amplitud; en el oscilador cuántico consiste en cambiar la ocupación de niveles separados por la energía del modo. Es útil separar, por tanto, sistema, variables, estado, modos disponibles y excitaciones efectivamente presentes.

![Figura E01.1. Un mismo lenguaje con objetos diferentes. En cada panel se identifican el sistema, el dominio, las variables, un modo y una excitación. El resorte puntual y el circuito concentrado no necesitan un campo espacial; la cavidad y el sólido sí admiten esa descripción. Son ejemplos ideales de elaboración propia.](figuras/E01_sistemas_campos.png){width=100%}

En el sólido, $\mathbf u(\mathbf r,t)$ describe desplazamientos de los átomos y un fonón es un cuanto de un modo vibratorio. Para describir un condensado superconductor se usa otro campo, $\Delta(\mathbf r,t)$, cuyo valor complejo tiene amplitud y fase. Una variación de su amplitud o de su fase es colectiva; una cuasipartícula electrónica es otro tipo de excitación. Qué variables se conservan depende de si queremos estudiar vibración, emparejamiento o respuesta eléctrica.

''' + e01[end:]
e01 = e01.replace('Los quarks y leptones son fermiones: un estado electrónico o fermiónico completamente especificado admite ocupación 0 o 1. Los bosones admiten otras ocupaciones.',
    'Los quarks y leptones son fermiones. Un **modo fermiónico completamente especificado**, incluido su espín, puede estar desocupado ($n=0$: ningún fermión en ese modo) u ocupado ($n=1$: un fermión en ese modo). Un orbital que admite dos orientaciones de espín contiene dos modos diferentes y puede alojar un electrón en cada uno. En un modo bosónico, como un modo de luz, las ocupaciones posibles son $n=0,1,2,\ldots$. El modo es la posibilidad disponible; su ocupación es parte de la descripción del estado.')
e01 = e01.replace('Masa, carga y espín no especifican por sí solos', 'La especie identifica propiedades compartidas: todos los electrones tienen la misma masa y carga. El espín es una propiedad de momento angular intrínseco. Masa, carga y espín no especifican por sí solos')
e01 = e01[:e01.index('**Conexión con A–D:**')] + '''**Fuentes:** [CERN, Modelo Estándar](https://home.cern/science/physics/standard-model/) y [MIT 8.03, modos normales y ondas](https://ocw.mit.edu/courses/8-03sc-physics-iii-vibrations-and-waves-fall-2016/). Los sistemas, el esquema y el ejercicio se construyen aquí para distinguir las descripciones.\n'''
(WORK/'E01.md').write_text(e01.strip()+'\n', encoding='utf-8')

e03 = section('# E.5.', '# E.6.')
e03 = e03.replace('y entender la envolvente de A.', 'y entender qué cambia al eliminar una variable interna estacionaria.')
e03 = e03.replace('Una función $F(y)$ recibe un número.', 'En el ejemplo más sencillo, una función $F(y)$ recibe un número.')
e03 = e03.replace('Se elige otro perfil', 'Variar la energía significa comparar cuánto vale **la misma regla de energía** para perfiles vecinos del mismo sistema. Se deforma su configuración; no se cambia arbitrariamente la fórmula ni se exige que esos perfiles sean una evolución temporal real.\n\nSe elige otro perfil')
e03 = e03.replace('El perfil ensayado es', 'Tomaremos $\\epsilon$ adimensional y $\\eta$ con las mismas unidades que $y$. El perfil ensayado es')
start = e03.index('$D\\mathcal F[y](\\eta)$ es un número')
end = e03.index('## E.5.2.')
e03 = e03[:start] + r'''La notación se lee en dos pasos: $D\mathcal F[y]$ es la regla de respuesta lineal calculada **en el perfil base** $y$; $(\eta)$ indica que esa regla se aplica **al perfil de perturbación completo** $\eta$. Los paréntesis admiten funciones como argumentos: no convierten a $\eta$ en una coordenada espacial. En $\eta(x)$ sí se evalúa el perfil en un punto $x$. Aquí se deriva respecto de $\epsilon$, manteniendo fijos ambos perfiles. El resultado $D\mathcal F[y](\eta)=J'(0)$ es un número, la pendiente de la energía en esa dirección.

Una versión con dos coordenadas aclara el papel de cada objeto. Si $U(\mathbf q)=(q_1^2+q_2^2)/2$, entonces

$$
DU[\mathbf q](\mathbf v)
=\left.\frac{d}{d\epsilon}U(\mathbf q+\epsilon\mathbf v)\right|_0
=q_1v_1+q_2v_2.
\tag{E03.a}
$$

Para $\mathbf q=(1,2)$ y $\mathbf v=(3,-1)$ la pendiente es $3-2=1$. $\mathbf v$ es el vector que elegimos para deformar la configuración, no el escalar respecto del cual derivamos. Al pasar de dos coordenadas a un perfil, la suma se convierte en una integral. Cuando, después de tratar los bordes, se puede escribir $D\mathcal F[y](\eta)=\int g(x)\eta(x)dx$, llamamos a $g(x)$ **derivada funcional**, escrita $\delta\mathcal F/\delta y(x)$. Esta función local y la pendiente integrada son objetos distintos.

''' + e03[end:]
e03 = e03.replace('Si los valores de $y$ en los extremos están fijados, las perturbaciones permitidas cumplen $\\eta(0)=\\eta(1)=0$. Con esa condición desaparece', r'''Fijar los extremos significa exigir los **mismos valores** $y_0$ y $y_1$ a todos los perfiles que comparamos. Esos valores pueden elegirse arbitrariamente al plantear el problema; una vez elegidos, no varían con $\epsilon$. Como $y(0)=y_0$ y $y(1)=y_1$, se exige

$$
\begin{aligned}
y_\epsilon(0)=y_0+\epsilon\eta(0)=y_0
&\ \Longrightarrow\ \epsilon\eta(0)=0,\\
y_\epsilon(1)=y_1+\epsilon\eta(1)=y_1
&\ \Longrightarrow\ \epsilon\eta(1)=0.
\end{aligned}\tag{E03.b}
$$

Estas igualdades deben valer para todo $\epsilon$ suficientemente pequeño, incluidos valores no nulos. Por eso $\eta(0)=\eta(1)=0$ aunque $y_0$ e $y_1$ no sean cero. Por ejemplo, una cuerda con extremos en alturas 2 y 5 puede ensayarse con $y(x)=2+3x$ y $\eta(x)=\sin(\pi x)$: la deformación cambia el interior y conserva ambas alturas. Con esa condición desaparece''')
e03 = e03.replace('En A hay variables internas que ya satisfacen su condición estacionaria.', 'Algunos modelos contienen variables internas que ya satisfacen su condición estacionaria.')
e03 = e03.replace('En A, resolver la variable espectral interna no obliga a que la amplitud del condensado esté también en equilibrio.', 'En un modelo con varias variables, resolver una variable interna no obliga a que las restantes estén también en equilibrio.')
e03 = e03.replace('**Conexión:** A aplica la envolvente; B distingue qué ocupaciones se fijan al derivar; C obtiene fuerzas y bordes de una energía espacial; D reúne esas ecuaciones. ', '')
(WORK/'E03.md').write_text(e03.strip()+'\n', encoding='utf-8')
print('Snapshot E-r01 and revised E01/E03 fragments ready.')
