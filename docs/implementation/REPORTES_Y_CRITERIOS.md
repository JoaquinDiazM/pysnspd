# Criterios prácticos de implementación y presentación

Aplicación desde la revisión del 24 de septiembre de 2026, por indicación del
usuario. El objetivo es incorporar el modelo SNSPD con evidencia útil para el
dispositivo, reutilizando métodos publicados y la infraestructura de la memoria.

## Selección del siguiente trabajo

Antes de añadir un algoritmo, comprobar si ya existe en producción, en el código
de los autores o en una biblioteca mantenida. Registrar la fuente y las unidades
que permiten reutilizarlo. Una publicación acredita el método que describe;
el acoplamiento concreto entre métodos requiere comprobar sus interfaces.

Cada cálculo debe resolver una decisión de implementación o una incertidumbre
de un observable previsto. Priorizar geometría dual, estados suaves resueltos,
bordes y puertos reales. Los casos límite sólo son necesarios cuando exponen un
riesgo efectivo de la implementación o del dominio físico que se va a usar.
La invariancia de unidades, orientación de flujos y ausencia de doble conteo
se verifican al conectar los bloques, sin reconstruir sus demostraciones.

Separar un fallo de ejecución de un criterio comparativo no satisfecho. Si ambas
trayectorias terminan, conservarlas y examinar qué señal incumple el criterio.
No convertir retrospectivamente el certificado fallido en pase ni pedir otra
campaña por defecto. Por indicación del usuario del 24 de septiembre, preferir
revisar una tolerancia innecesariamente exigente y cerrar la discrepancia con
una nueva decisión de aceptación práctica, en lugar de aparcarla. Conservar
intacto el resultado bajo el criterio anterior y explicar la escala y finalidad
que justifican el margen nuevo. Esto no acredita más precisión que la medida.
Antes de la próxima corrida declarar
qué observables deciden la aceptación y cuáles son diagnósticos secundarios.

Una comparación temporal sólo aporta refinamiento independiente en los tiempos
donde las secuencias de pasos difieren. Un porcentaje referido a una señal casi
nula se acompaña de su diferencia absoluta y de una escala física pertinente
del mismo observable. Ese contexto no reemplaza el resultado original.

La revisión no requiere esperar a una discrepancia catastrófica para reparar
un error real: unidades incorrectas, estados no finitos, crecimiento inestable,
violaciones materiales de conservación o términos físicos ausentes no se
corrigen aumentando un umbral. Si el desacuerdo es pequeño frente al objetivo,
registrar su aceptación al margen revisado y continuar sin otra campaña por
defecto. El 2 % adoptado para los controles débiles de esta iteración es un
presupuesto de desarrollo; no es una precisión universal de latencia ni un
error experimental del material.

Se prioriza probar el paso KWT heredado. Su ecuación local cuadrática determina
la amplitud del paso nuevo, pero el avance sigue siendo Euler de primer orden.
Comparar los observables relevantes a paso y medio paso permite decidir su
utilidad, sin rediseñar el integrador por anticipado.

## Qué debe decir cada figura

Cada título o pie identifica:

1. El caso físico: control térmico, ensayo del circuito o detector; geometría,
   temperatura, excitación y tiempo observados.
2. El campo o magnitud representada, su fórmula o procedimiento de extracción
   y sus unidades. Una norma de campo complejo no se llama amplitud puntual.
3. Si es campo total, diferencia respecto a una referencia o diferencia entre
   dos resoluciones. Identificar ambas soluciones y la referencia sustraída.
4. El denominador de cualquier normalización. Para una norma espacial,
   declarar pesos de área o de arista y la región integrada.
5. Qué significan colores y trazos. Mapas comparables usan la misma escala;
   curvas solapadas se acompañan de un panel de diferencias cuando corresponda.

Por ejemplo, `||d_A(t)-d_b(t)||_M/||d_A(0)-d_b(0)||_M`, con
`d=Delta/(k_B Tc)` y `||z||_M²=sum_i m_i |z_i|²`, es una respuesta espacial
normalizada y adimensional. No es `|Delta|` de un punto ni un voltaje.
`Im(conj(d_i)G_i)/m_i` es una fuerza variacional de fase en las unidades de la
acción; no es un torque mecánico. Los gráficos deben explicitar esa distinción.

Los datos guardados se muestran como muestras. Una línea que los une guía la
lectura y no constituye una nueva solución interpolada. El informe no presenta
un control sin fotón como hotbelt, jitter, latencia o V_out del detector.
