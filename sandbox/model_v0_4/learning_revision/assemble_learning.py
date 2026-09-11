"""Assemble independent reviewed lesson fragments into the stable 0.4 notebook."""
from pathlib import Path
import json
import re
import hashlib

ROOT=Path(__file__).resolve().parents[3]
DOCS=ROOT/'docs/modelo_v0_4'
WORK=Path(__file__).resolve().parent
ORDER=['E01','E04','E06','E07','E03','E05','E09','E10']
intro=r'''---
title: "E. Cuaderno de aprendizaje del modelo"
subtitle: "Clases autónomas · Edición del modelo 0.4"
date: "9 de septiembre de 2026 · Cuaderno E-r02"
lang: es
---

\setlength{\parskip}{4pt plus 1pt}

# E.0. Ruta y registro de aprendizaje

Cada clase parte de un sistema concreto, presenta sus variables y desarrolla un ejemplo antes de proponer actividades. Puede leerse de manera aislada: los conceptos necesarios se recuperan en el propio texto. La tabla ofrece material previo recomendado, sin darlo por dominado ni exigir una lectura lineal. Las relaciones entre clases se expresan al volver a usar las ideas; los códigos se reservan para orientarse y registrar respuestas.

| ID estable | Clase | Estado | Material previo recomendado |
|:--|:--|:--|:--|
| E01 | Sistemas, campos, estados y excitaciones | ABIERTO | Posición y velocidad; vectores; corriente y tensión |
| E04 | Electrones y huecos: materiales y fronteras | ABIERTO | E01; carga eléctrica; corriente convencional; energía potencial |
| E06 | Ocupaciones y operadores | ABIERTO | E01, E04; matrices; derivadas; mecánica de una masa y un resorte |
| E07 | Una columna de Nambu | ABIERTO | E04, E06; números complejos; valores propios |
| E03 | Variar una energía funcional | ABIERTO | E01; derivadas; regla del producto y de la cadena; integración por partes |
| E05 | Modos, distribuciones y DOS | ABIERTO | E01; oscilador armónico; integrales; valores propios; probabilidad básica |
| E09 | DFT: densidad y estructura electrónica | ABIERTO | E01, E03, E04, E05; normalización; potencial y energía; orbitales |
| E10 | DFPT: respuesta e interacción electrón–fonón | ABIERTO | E03, E05, E09; E07 para el factor superconductor; Taylor; unidades y conservación de energía |

Los números E01, E03, etc. son identificadores permanentes; E.1, E.2, etc. sólo indican el orden de lectura. E02 está descartada y fue sustituida por E06/E07. E08 se retira y se sustituye por las nuevas E09 (DFT) y E10 (DFPT). No se asignan actividades ni se recomienda estudiar las clases retiradas. Su retiro es una decisión editorial, no una evaluación del aprendizaje.

No hay respuestas entregadas ni calificaciones asignadas. Cada clase propone tres actividades, por 4, 3 y 3 puntos; se cierra con al menos 8 puntos y su condición esencial satisfecha. Los ejemplos resueltos son distintos de esas actividades. Las observaciones de lectura sirven para revisar la enseñanza y no se califican como respuestas incorrectas. Las respuestas futuras y su evidencia se conservan aunque se revise una clase; si cambia una consigna, se conserva también su versión anterior.

**Revisión de esta entrega:** E01, E03, E04, E05, E06 y E07 pasan de r1 a r2; E09 y E10 comienzan en r1. El cuaderno ensamblado es E-r02 y sigue perteneciendo a la edición del modelo 0.4. Una corrección posterior incrementará sólo la revisión de las clases afectadas y la del ensamblado; no cambiará automáticamente A–D ni la versión del modelo. Las ecuaciones antiguas conservan sus etiquetas; las nuevas usan el código de su clase para que añadir contenido no renumere las demás. El historial y la regla completa se guardan junto al documento editable.
'''

parts=[intro.strip()]
for code in ORDER:
    text=(WORK/f'{code}.md').read_text(encoding='utf-8').strip()
    text=re.sub(r'\s*\\Needspace\{12\\baselineskip\}\s*$', '', text)
    text=re.sub(r'^\\Needspace\{12\\baselineskip\}\s*', '', text)
    rev=1 if code in ['E09','E10'] else 2
    text=text.replace('**Estado: ABIERTO.',f'**Revisión: r{rev}. Estado: ABIERTO.',1)
    separator='\\clearpage' if code=='E01' else '\\Needspace{14\\baselineskip}\n\n\\bigskip'
    parts.append(separator+'\n\n'+text)
parts.append(r'''\Needspace{14\baselineskip}

\bigskip

# E.9. Respuestas y próxima revisión

Hay 24 respuestas pendientes, tres por cada una de las ocho clases activas. No se requiere resolver actividades de E02 o E08. Las respuestas se escriben bajo su identificador en el documento editable; la pauta y las calificaciones siguen pendientes de una entrega.

**Comentarios de aprendizaje:** _Indicar qué símbolo, concepto o paso necesita una explicación adicional._

La siguiente revisión conservará las respuestas y modificará sólo el material que las observaciones justifiquen. El estado de aprendizaje y la revisión del contenido se registran por separado: una clase revisada no se considera automáticamente aprobada.
''')
source='\n\n'.join(parts)+'\n'
assert source.count('**Respuesta ')==24
assert len(re.findall(r'^# E\.',source,re.M))==10
tags=re.findall(r'\\tag\{([^}]+)\}',source)
assert len(tags)==len(set(tags))
target=DOCS/'E_cuaderno_de_aprendizaje_v0_4.md'
stamp=WORK/'assembly_source.json'
if target.exists() and target.read_text(encoding='utf-8') != source:
    previous=json.loads(stamp.read_text(encoding='utf-8')) if stamp.exists() else {}
    actual=hashlib.sha256(target.read_bytes()).hexdigest()
    if previous.get('sha256') != actual:
        raise RuntimeError('El Markdown contiene una edición independiente. Es la fuente canónica: no se sobrescribe con fragmentos históricos. Usar build_learning.py directamente.')
target.write_text(source,encoding='utf-8')
stamp.write_text(json.dumps({'notebook_revision':'E-r02','sha256':hashlib.sha256(target.read_bytes()).hexdigest(),
    'purpose':'Guarda de montaje histórico. Editar el Markdown canónico y compilar con build_learning.py.'},indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
records={'model_edition':'0.4','notebook_revision':'E-r02','date':'2026-09-09',
 'classes':{code:{'revision':1 if code in ['E09','E10'] else 2,'learning_state':'ABIERTO'} for code in ORDER},
 'retired':{'E02':{'revision':1,'replaced_by':['E06','E07']},'E08':{'revision':1,'replaced_by':['E09','E10']}},
 'activities_pending':24,'previous_snapshot':'historial/E-r01',
 'notes':'Las revisiones son editoriales. No certifican aprendizaje ni validación física del detector.'}
(DOCS/'E_registro_de_revisiones.json').write_text(json.dumps(records,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(f'Assembled {len(ORDER)} lessons, {len(tags)} equations, 24 unanswered activities.')
