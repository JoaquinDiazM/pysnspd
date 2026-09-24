# Etapa 4: cierre de desarrollo con frontera física explícita

La campaña de dos referencias con corriente y cinco respuestas acopladas terminó.
Se cierra el **desarrollo sin fotón de etapa 4**. No se declara completado todo
el contrato dinámico ni se admite todavía el modelo actualizado para producción.

El [informe final](Informe_cierre_etapa_4.md)
([PDF](../../../../output/pdf/implementation/Informe_cierre_etapa_4.pdf)) muestra
campos y observables físicos: supresión del gap por corriente, relajación de
amplitud y fase, respuesta espacial, inductancia y señal del circuito completo.

## Resultado que permite continuar

- Dos estados superconductores a 0,9 K, sobre 160 x 80 nm y 1712 nodos:
  3,368 y 13,002 microamperios. El mínimo del gap pasa de 0,99849 a 0,97471
  de la referencia sin corriente; no aparece una región normal.
- La respuesta compleja refinada cambia 1,589 % en admitancia y 0,0005825 %
  en Vout, con el circuito completo y la misma excitación de fuente.
- El control temporal previo de 1 ps conserva la malla dual y Euler KWT
  heredado. No es necesario repetirlo ni cambiar el integrador.

## Frontera que el cierre no elimina

La respuesta es predominantemente inductiva. Su estabilidad no acredita el
calentamiento: la parte real de la admitancia cambia 45,66 % respecto al
refinado, aunque su diferencia absoluta sólo sea 0,8194 % de la admitancia
total. El calor radial candidato supera en 31,37 % la potencia de puerto en
el caso refinado. Esta comparación no incluye un balance independiente de
trabajo DC, reservorios y energía interna; no se convierte en calor validado.

La corriente AC integrada a través de cortes internos se aparta hasta 4,045 %
de la corriente del puerto izquierdo en el refinado. La coincidencia de los
dos terminales no basta para garantizar continuidad en todo el interior.
Se admite la aproximación reactiva para desarrollo, sin asignar 2 % de precisión
a la disipación ni a la conservación local de todos los casos.

La [decisión estructurada](closure_decision.json) conserva estas distinciones.
La [entrada a etapa 5](../../stage5/README.md) comienza por completar las
interfaces que afectan a calor y corriente antes de un transiente fotónico.
No hay otra corrida larga solicitada para cerrar esta entrega.

## Reproducibilidad y almacenamiento

`data/` contiene unos 3 MB de campos finales y recibos revisados. Las referencias
NPZ son extractos para análisis, **no puntos de reinicio del solver**. Los campos Matsubara de las referencias, momentos integrados y salidas
originales permanecen en
`/home/jdiaz/scratch/stage4_final_coupling_20260924`; su procedencia y los hashes
se registran en [data_manifest.json](data_manifest.json). Se conserva el plan
ejecutado, el código y las fuentes originales de los cálculos anteriores.

Desde la raíz del repositorio, con el entorno científico activo:

```bash
python sandbox/stage4_core/analyze_final_closure.py
python sandbox/stage4_core/plot_final_closure.py
python sandbox/stage4_core/build_final_closure_report.py
```

Estos comandos sólo posprocesan datos. La [limpieza del repositorio](../../ARCHIVOS_ARCHIVADOS.md)
retira copias exactas e informes intermedios, manteniendo cierres, fuentes,
datos canónicos y restauración verificable. No reescribe Git ni modifica v1.0.0.
