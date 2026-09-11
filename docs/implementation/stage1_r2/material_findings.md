La comparación es de **DOS fonónica**, no de la DOS electrónica de Usadel. Se invocó el loader real del legado con los argumentos usados por PRE.

El legado conserva la parte positiva y recorta 361 valores negativos. La integral pasa de 0.708420058 a 0.719950139; no hay renormalización. Conserva cuatro pares de abscisas repetidas; el ordenamiento invierte un par y modifica la integral bruta en 1.60e-10, antes del recorte. α²F coincide punto a punto después del mismo ordenamiento: por eso sus dos trazos se superponen. La nueva figura muestra las coincidencias mediante estilos distintos y separa el residual y la cola para visualizar el cambio.

**R2 mantiene la entrada rechazada.** No existe una nueva curva fonónica experimental que superponer. La procedencia numérica se verifica contra [el archivo público fijado por revisión](https://github.com/qnngroup/proj-KE-solver/blob/5b6bd747f80016da5ccd51db73c110a8ecc6abf6/nbn-a2f-ph.dat), que carece de la cabecera añadida localmente. La unidad y la base de normalización de la ordenada siguen sin certificar.

El [solver de la fuente](https://github.com/qnngroup/proj-KE-solver/blob/5b6bd747f80016da5ccd51db73c110a8ecc6abf6/solver.m#L93-L113) referencia `nbn-a2f-ph_2.dat`, ausente del árbol publicado consultado. Su conversión de unidades no demuestra que este archivo sin sufijo tenga las mismas convenciones. No se adopta una reinterpretación por conveniencia.

La verificación independiente R2 pasa 25 controles y reproduce la integral de la primera etapa con diferencia 1.11e-16. Tiempo de ejecución: 1.222 s.
