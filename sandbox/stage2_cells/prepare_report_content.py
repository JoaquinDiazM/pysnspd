"""Assemble the concise stage-2 report from adjudicated evidence, never pilots."""
from pathlib import Path
import argparse
import json
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'docs/implementation/stage2'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def text(value, kind='body'):
    return dict(type=kind,text=value)


def figure(name,caption):
    return dict(type='image',path=f'docs/implementation/stage2/figures/{name}.png',caption=caption)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--one', type=Path, required=True)
    parser.add_argument('--two', type=Path, required=True)
    parser.add_argument('--verdict-text', required=True)
    parser.add_argument('--tests-text', required=True)
    args = parser.parse_args()
    one, two = read(args.one), read(args.two)
    a1, a2 = one['final']['amplitudes'], two['final']['amplitudes']
    ledger=max(one['energy_ledger_scaled_max'],two['energy_ledger_scaled_max'])
    instant=max(one['instantaneous_residual_max'],two['instantaneous_residual_max'])
    nodes=two['occupation_mesh']['electron_states']
    pages=[dict(title='Etapa 2 | Celdas acopladas',blocks=[
        text(args.verdict_text),
        text('Resultado principal','subheading'),
        text('El condensado, las ocupaciones electrónicas y las fonónicas evolucionan juntos. '
             'Las mismas reacciones intercambian energía entre electrones y fonones; '
             'el transporte conecta celdas con espectros distintos. La disipación del '
             'condensado, la fuente y el escape se contabilizan una sola vez.'),
        dict(type='table',widths=[.39,.36,.25],rows=[['Ensayo','Resultado final','Interpretación'],
            ['Una celda aislada',f'Amplitud: 0,6000 → {a1[0]:.4f}','Recuperación con cinética activa'],
            ['Dos celdas',f'Amplitudes: {a2[0]:.4f} / {a2[1]:.4f}','Intercambio entre espectros distintos'],
            ['Balance energético',f'Error escalado ≤ {ledger:.2e}','Límite: 1e-7'],
            ['Balance instantáneo',f'Residuo ≤ {instant:.2e}','Límite: 1e-10'],
            ['Resolución electrónica',str(nodes)+' estados por celda','Malla complementaria explícita']]),
        text(args.tests_text),
        text('Alcance de los resultados','subheading'),
        text('Los tiempos y materiales de estos ensayos son sintéticos. Se ha comprobado '
             'el sistema reducido de celdas; las tasas absolutas de NbN, los gradientes '
             'espaciales del condensado, los reservorios, Poisson y el circuito requieren '
             'sus propias pruebas. El solver de producción y la versión v1.0.0 se conservan.'),
    ]),dict(title='Evolución simultánea y balance',blocks=[
        figure('coupled_dynamics','Mismo sistema de ecuaciones para las curvas y los registros de energía. '
               'Una transferencia electrón-fonón negativa significa energía neta de fonones a electrones. '
               'La curva de transporte representa intercambio interno y no se suma otra vez al balance global.'),
        text(f"En dos celdas, el transporte acumulado es {two['final']['transport']:.5f} y el intercambio "
             f"neto electrón-fonón es {sum(two['final']['electron_to_phonon']):.5f}, en unidades de N₀Δ₀². "
             'Los canales están activos durante la recuperación de la amplitud.'),
        text('La temperatura equivalente es una lectura de la energía electrónica de cada celda. '
             'El sistema integra las ocupaciones, por lo que puede conservar una distribución no térmica.'),
        text('Escalas: energía de partícula en Δ₀, densidad de energía en N₀Δ₀² y tiempo en '
             't_ref = 1 ps de ensayo. El residuo del balance se escala por el máximo de 1, '
             'la energía inicial absoluta y la energía externa inyectada.','small'),
    ]),dict(title='Qué cambia en las distribuciones',blocks=[
        figure('resolved_populations','Ensayo aislado. Cada distribución electrónica se dibuja en su propio '
               'espectro; la diferencia inferior compara la población final con Fermi-Dirac a igual energía. '
               'Los paneles de diferencias hacen visibles cambios que quedarían ocultos por curvas superpuestas.'),
        text('El calentamiento y la relajación no obligan a reemplazar p por una distribución térmica. '
             'Las reacciones y BGK determinan cuánto se aproxima a ella. Las ocupaciones fonónicas '
             'se redistribuyen simultáneamente, conservando el intercambio energético compartido.'),
    ]),dict(title='La resolución se decide con distribuciones',blocks=[
        figure('transport_convergence','Comparación con una referencia continua independiente. '
               'El panel de residuos distingue curvas casi superpuestas; el error local del perfil '
               'no debe confundirse con la norma integrada usada en el criterio.'),
        text('La malla original de 180 estados conservaba energía, pero no resolvía suficientemente '
             'ciertas distribuciones. Se adoptó una malla complementaria y una partición de eventos '
             'que sigue tanto los nodos electrónicos como las energías fonónicas. '
             'Los pilotos fallidos permanecen registrados.'),
        text('La región cercana al borde espectral se verifica también mediante el momento conjugado '
             'a la corriente. Su resolución es necesaria para la siguiente etapa, aunque las '
             'trayectorias principales de esta entrega tengan Gamma = 0.'),
    ]),dict(title='Verificación y decisión de continuidad',blocks=[
        text('La evidencia está organizada en cinco controles independientes.','body'),
        dict(type='table',widths=[.31,.69],rows=[['Control','Qué se exige'],
            ['Identidades discretas','Cada evento conserva energía; dispersión y transporte conservan número de cuasipartículas.'],
            ['Equilibrio y fronteras','Cancelación Fermi/Bose, 0 ≤ p ≤ 1, n ≥ 0, sin recortar ocupaciones.'],
            ['Consistencia continua','Potencias y respuestas poblacionales separadas, con referencias independientes y refinamientos de malla y cuadratura.'],
            ['Integración temporal','Tres pasos temporales y un integrador de referencia; balance y errores intermedios además del estado final.'],
            ['Soporte y fuerzas','Cortes infrarrojo/superior, causalidad y derivadas de la misma energía en los campos visitados.']]),
        text('Los controles negativos detectan el doble conteo de calor y la omisión del trabajo '
             'espectral. Un error de normalización de las tasas también es rechazado. '
             'Conservar energía, por sí solo, no basta para superar estos controles.'),
        text('Siguiente etapa','subheading'),
        text('Seguir el contrato NEXT_STAGE.md: comprobar primero la energía y rigidez espacial, '
             'después interfaces y bordes, y luego Poisson y circuito. La base de celdas '
             'reduce la incertidumbre cinética, pero no sustituye la admisión de esos términos.'),
        text('Trazabilidad','subheading'),
        text('Dictamen: stage2_admission.json. Criterios: acceptance_criteria.json y anexo de malla '
             'complementaria. Detalle de iteraciones: NUMERICAL_DECISIONS.md. '
             'Fuentes del modelo: documentos B y D de la versión 0.4; espectro y energía común, A y C. '
             'La carpeta review contiene las referencias independientes y sus errores. '
             'Los comandos reproducibles se conservan en GEMINGA_COMMANDS.md.','small'),
    ])]
    content=dict(title='pySNSPD | Resultados de la etapa 2',
        subtitle='Celdas superconductoras con cinética y energía compartida · Modelo 0.4 · 21 septiembre 2026',pages=pages)
    (DATA/'report_content.json').write_text(json.dumps(content,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')


if __name__ == '__main__':
    main()
