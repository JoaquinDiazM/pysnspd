# pySNSPD template repository

Plantilla modular para reconstruir la librería de simulación multiescala de SNSPDs de forma ordenada, verificable y extensible.

Esta versión **no implementa todavía la física numérica**. Todos los módulos contienen funciones con `return 0` y docstrings extensos para fijar responsabilidades, interfaces esperadas y orden de trabajo. La intención es que cada bloque pueda desarrollarse y probarse por separado antes de volver a acoplar todo el pipeline.

---

## 1. Filosofía física del framework

El modelo separa el problema en cuatro niveles:

1. **Material y espectro microscópico**  
   Se obtiene una descripción de superconductor sucio mediante Usadel. Este bloque entrega cantidades como
   `Delta`, `rho(E; |Delta|, q)`, `D`, `sigma_n`, corriente crítica estimada, relación corriente--momento condensado y catálogos espectrales.

2. **Proyección cinética QP--fonón**  
   A partir de los kernels tipo Simon/MIT se proyectan las ecuaciones cinéticas sobre energía. El resultado no es resolver siempre `f(E,t)` y `n(Omega,t)` completos, sino construir potencias microscópicas:

   - `P_ep_S`: potencia de scattering QP--fonón.
   - `P_ep_R`: potencia de recombinación / pair-breaking.
   - `P_esc`: escape fonónico hacia sustrato.

   Estas potencias se tabulan mediante integrales de phase space `mathcal{J}_S` y `mathcal{J}_R`, dependientes de `Te`, `Tph`, `|Delta|` y `q`.

3. **Dinámica mesoscópica gTDGL**  
   La dinámica del parámetro de orden `Psi = R exp(i phi)` decide si la perturbación térmica se convierte en un evento disipativo. Este bloque resuelve amplitud, fase, potencial eléctrico, corrientes y continuidad.

4. **Circuito externo y lectura**  
   La caída de voltaje interna `V_TDGL(t)` se acopla a un circuito mínimo o extendido para obtener `I_SNSPD(t)` y `V_out(t)`.

La contribución conceptual del framework es que no usa una ecuación térmica `2T` como receta aislada. La ecuación térmica se interpreta como una **proyección energética** de la cinética QP--fonón, compatible con gTDGL, Usadel y el circuito.

---

## 2. Etapas obligatorias de una simulación

El repositorio queda organizado alrededor de tres etapas principales.

### 2.1 PRE-run

La etapa PRE-run construye todos los objetos caros, reutilizables y casi-estáticos.

Debe generar:

- malla bidimensional de la nanocinta;
- triangulación tipo Delaunay;
- lista de aristas y conectividad;
- geometría de contactos y bordes;
- parámetros materiales derivados;
- solución Usadel para el estado de bias;
- catálogo fino de DOS de cuasipartículas `rho(E; |Delta|, q)`;
- catálogo fino de integrales de phase space `mathcal{J}_S` y `mathcal{J}_R`;
- metadatos suficientes para reproducibilidad.

La PRE-run debe poder paralelizarse, porque la construcción de catálogos finos puede ser costosa. La configuración debe permitir decidir cuántos workers, procesos o nodos se usarán.

### 2.2 SS-run

La etapa SS-run busca la condición inicial estacionaria previa al fotón.

Debe cargar los resultados de PRE-run y construir una predicción analítica inicial razonable:

- `Te = Tph = T_bath`;
- amplitud `|Delta|` cercana al valor Usadel/BCS compatible con `T_bias` e `I_bias`;
- fase con rampa longitudinal consistente con el momento condensado `q`;
- potencial eléctrico inicialmente cercano a cero;
- corriente total compatible con `I_bias`;
- circuito en estado estacionario.

El objetivo no es partir de ruido ni de una condición arbitraria, sino desde una aproximación física que haga converger rápido hacia un estacionario numéricamente limpio.

La salida de SS-run debe ser la condición inicial oficial de PHOTON-run.

### 2.3 PHOTON-run

La etapa PHOTON-run evoluciona el sistema después de la absorción del fotón.

Debe cargar:

- el estado estacionario de SS-run;
- catálogos de PRE-run;
- configuración del fotón;
- configuración térmica y de circuito.

La perturbación inicial se introduce como una phonon bubble construida desde la energía del fotón. Luego se evoluciona:

- temperatura electrónica efectiva `Te(x,y,t)`;
- temperatura fonónica efectiva `Tph(x,y,t)`;
- parámetro de orden `Psi(x,y,t)`;
- potencial eléctrico;
- corrientes;
- circuito externo;
- señales `V_TDGL(t)`, `I_SNSPD(t)`, `V_out(t)`.

---

## 3. Manejo obligatorio de archivos

Esta plantilla fuerza una idea central: **los datos crudos no deben vivir dentro de la librería**.

El usuario debe configurar explícitamente un directorio externo:

```yaml
project:
  big_data_root: /path/to/my_big_data
```

Ese directorio se usará como base para toda acumulación de datos pesados:

```text
big_data_root/
  raw/
    <run_name>/
      manifest.yaml
      pre/
      ss/
      photon/
  plots/
    <run_name>/
      manifest.yaml
      figures/
  logs/
    <run_name>/
```

El nombre de la run debe ser idéntico en datos crudos y plots. Por ejemplo:

```text
raw/run_035uA_1064nm/
plots/run_035uA_1064nm/
logs/run_035uA_1064nm/
```

Por diseño, **es mala idea renombrar una run después de simular datos o generar plots**, porque el módulo de plotting debe inferir la relación entre datos y figuras a partir del mismo `run_name`.

La librería debe guardar siempre:

- configuración completa usada en la run;
- versión del paquete;
- fecha y hora;
- hash o resumen del input material cuando corresponda;
- parámetros de malla;
- parámetros de catálogos;
- parámetros térmicos;
- parámetros gTDGL;
- parámetros de circuito;
- rutas absolutas relevantes.

---

## 4. Estructura del repositorio

```text
pysnspd_template_repo/
  README.md
  pyproject.toml
  .gitignore
  configs/
    example_project.yaml
  pipelines/
    00_configure_project.py
    01_prerun_template.py
    02_ss_run_template.py
    03_photon_run_template.py
    04_plot_run_template.py
  pysnspd/
    __init__.py
    config.py
    io/
      manager.py
    mesh/
      delaunay.py
      edges.py
    usadel/
      parameters.py
      solver.py
      catalog.py
    kinetic/
      eliashberg.py
      phase_space.py
      powers.py
      thermal.py
      phonon_bubble.py
    gtdgl/
      initial_guess.py
      solver.py
      currents.py
      poisson.py
      boundary.py
    circuit/
      bias_tee.py
    runs/
      prerun.py
      ssrun.py
      photonrun.py
      orchestrator.py
    plotting/
      load.py
      figures.py
    validation/
      checks.py
    utils/
      parallel.py
      logging.py
      units.py
  tests/
    test_imports.py
```

---

## 5. Responsabilidad de cada módulo

### `pysnspd.config`

Define la interfaz de configuración global. Más adelante debe validar campos obligatorios, rutas, unidades y coherencia física.

### `pysnspd.io.manager`

Debe centralizar el manejo de archivos. Ningún solver debería construir rutas a mano. Todo acceso a datos crudos, plots, logs y manifests debe pasar por este módulo.

### `pysnspd.mesh`

Contiene malla, Delaunay, aristas, conectividad, bordes y geometría de contactos.

### `pysnspd.usadel`

Contiene calibración material, ecuación Usadel, barridos en `q`, relación corriente--momento, DOS de QPs y catálogos espectrales.

### `pysnspd.kinetic`

Contiene funciones relacionadas con Eliashberg, phase-space integrals, potencias proyectadas, cierre térmico y phonon bubble.

### `pysnspd.gtdgl`

Contiene predicción analítica inicial, solver gTDGL, corrientes, Poisson y condiciones de borde.

### `pysnspd.circuit`

Contiene el circuito externo, inicialmente bias tee / carga de 50 Ohm / ecuaciones para `I_SNSPD` y `V_out`.

### `pysnspd.runs`

Contiene wrappers de alto nivel para PRE-run, SS-run y PHOTON-run. Estos módulos no deben tener detalles matemáticos internos: deben orquestar los módulos inferiores.

### `pysnspd.plotting`

Debe cargar datos desde `big_data_root/raw/<run_name>` y guardar figuras en `big_data_root/plots/<run_name>`. No debe modificar datos crudos.

### `pysnspd.validation`

Debe reunir chequeos de conservación de corriente, consistencia energética, ausencia de voltaje espurio, continuidad de fase, estabilidad del estacionario y compatibilidad entre catálogos.

---

## 6. Roadmap metodológico de implementación

El desarrollo de `pySNSPD` se organiza como una secuencia de objetivos verificables. La idea no es implementar todo el modelo de una vez, sino construir una cadena modular donde cada bloque pueda probarse de forma independiente antes de acoplarse al siguiente. La marca ★ indica el objetivo actualmente en desarrollo.

### Objetivos generales

**OG1. Infraestructura reproducible de simulación.**
Construir una base de librería que permita configurar proyectos, organizar datos pesados fuera del repositorio, registrar metadatos de cada corrida y reproducir resultados sin depender de rutas o nombres ambiguos.

**OG2. Construcción de catálogos microscópicos y mesoscópicos.**
Precomputar los objetos costosos del modelo, incluyendo malla, aristas, parámetros Usadel, densidad de estados de cuasipartículas y catálogos de integrales de fase espacial (\mathcal J_S) y (\mathcal J_R).

**OG3. Simulación acoplada de detección.**
Integrar condición estacionaria, excitación por fotón, evolución térmica, dinámica gTDGL, conservación de corriente y circuito externo para estudiar cuándo una perturbación local produce un evento disipativo observable.

### Objetivos específicos y resultados esperados

1. **Configurar proyecto y manejo de archivos. ★**
   Resultado esperado: lectura de configuración, validación de `big_data_root`, creación automática de carpetas `raw/`, `plots/`, `logs/`, `catalogs/` y escritura de manifests por corrida.

2. **Implementar malla, aristas, bordes y visualización básica.**
   Resultado esperado: generación reproducible de mallas Delaunay, identificación de bordes/contactos, aristas internas y plots diagnósticos de geometría.

3. **Implementar bloque Usadel y catálogo DOS.**
   Resultado esperado: catálogo (\rho(E;|\Delta|,q)), relación corriente--momento superconductivo y parámetros materiales derivados para la etapa mesoscópica.

4. **Implementar catálogos (\mathcal J_S) y (\mathcal J_R).**
   Resultado esperado: tablas interpolables de integrales de scattering y recombinación/pair-breaking compatibles con la proyección energética del Apéndice A.

5. **Construir el handler de PRE-run.**
   Resultado esperado: etapa reproducible y paralelizable que genere malla, catálogos microscópicos, metadatos y rutas de salida antes de cualquier evolución temporal.

6. **Construir condición inicial analítica para SS-run.**
   Resultado esperado: campo inicial razonable para (|\Delta|), fase, corriente y potencial, reduciendo el tiempo necesario para alcanzar un estado estacionario numérico.

7. **Implementar gTDGL estacionario sin fotón.**
   Resultado esperado: solver estacionario capaz de relajar la nanocinta polarizada sin generar artefactos eléctricos dominantes.

8. **Validar conservación de corriente.**
   Resultado esperado: diagnósticos de (\nabla\cdot\mathbf j), corriente integrada en contactos, fase desenvuelta y voltaje espurio en estado estacionario.

9. **Implementar phonon bubble.**
   Resultado esperado: fuente inicial dependiente de la energía del fotón y de la geometría local de absorción.

10. **Implementar PHOTON-run térmico desacoplado.**
    Resultado esperado: evolución de (T_e) y (T_{ph}) usando (P_{ep}^{S}), (P_{ep}^{R}), difusión y escape, sin acoplar todavía la respuesta gTDGL.

11. **Acoplar PHOTON-run térmico con gTDGL.**
    Resultado esperado: evolución conjunta de temperatura electrónica, parámetro de orden, corriente y potencial después de la absorción.

12. **Acoplar circuito externo.**
    Resultado esperado: evolución de (I_{\rm SNSPD}(t)), (V_{\rm TDGL}(t)) y (V_{\rm out}(t)) con una lectura compatible con una línea de (50,\Omega).

13. **Implementar comparación entre corridas.**
    Resultado esperado: módulo de ploteo que use el mismo `run_name` para datos y figuras, permitiendo comparar bias, energía del fotón, geometría, catálogos y respuesta eléctrica.


---

## 7. Convención de nombres

Recomendación para `run_name`:

```text
<material>_<width>nm_<Ibias>uA_<photon>nm_<tag>
```

Ejemplo:

```text
NbN_120nm_35uA_1064nm_test001
```

El mismo nombre debe usarse para:

```text
raw/NbN_120nm_35uA_1064nm_test001/
plots/NbN_120nm_35uA_1064nm_test001/
logs/NbN_120nm_35uA_1064nm_test001/
```

---

## 8. Estado actual de esta plantilla

Esta plantilla solo fija arquitectura. Las funciones existen, pero todavía no calculan física real. La mayoría tiene esta forma:

```python
def some_function(...):
    """Describe qué hará la función cuando sea implementada."""
    return 0
```

Esto permite empezar de manera limpia: primero se estabiliza el repositorio, después se llenan los bloques uno por uno.

---

## 9. Ejecución esperada futura

La secuencia futura será:

```bash
python pipelines/00_configure_project.py --config configs/example_project.yaml
python pipelines/01_prerun_template.py --config configs/example_project.yaml --run-name NbN_120nm_35uA_1064nm_test001
python pipelines/02_ss_run_template.py --config configs/example_project.yaml --run-name NbN_120nm_35uA_1064nm_test001
python pipelines/03_photon_run_template.py --config configs/example_project.yaml --run-name NbN_120nm_35uA_1064nm_test001
python pipelines/04_plot_run_template.py --config configs/example_project.yaml --run-name NbN_120nm_35uA_1064nm_test001
```

Por ahora, estos scripts solo llaman funciones placeholder.

---

## 10. Principio de diseño

Cada etapa debe ser reproducible, cacheable y testeable. Una PHOTON-run no debería recalcular una DOS si ya existe una PRE-run compatible. Una SS-run no debería depender de archivos sueltos no registrados. Un script de plotting no debería modificar datos crudos. El `run_name` y el `manifest` son el pegamento entre todo.
