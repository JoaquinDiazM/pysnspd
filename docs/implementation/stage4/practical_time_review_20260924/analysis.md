# Diagnostico de la corrida termica no lineal

Las dos trayectorias llegaron a 1 ps en 1089.49 s. La interrupcion fue el control posterior de refinamiento: 95/96 comparaciones cumplen y 1 no. No fallo Newton ni la integracion.
Hasta 0.1 ps ambas salidas son identicas porque utilizaron los mismos pasos. De las 24 comparaciones en tiempos con trayectorias diferentes (0.3 y 1 ps), 23 cumplen. No son 96 pruebas independientes ni un estudio de orden de convergencia.

## Comparacion que no cumple

- Campo: torque variacional de fase inducido por la perturbacion radial de amplitud, restando el nucleo base al mismo tiempo; t=1 ps.
- Diferencia primaria-refinada: 8.3046528e-08; tolerancia registrada: 3.76818628e-08.
- Escala inicial de ese torque: 5.53637256e-06; diferencia: 1.50002% de esa escala.
- Torque restante refinado: 1.76208015e-08; diferencia: 471.298% de la senal restante.
- Diferencia contextual respecto al torque inicial principal (sonda angular): 0.00456644%. Esta escala adicional no cambia el resultado registrado.

## Alcance fisico de la evidencia

El torque de fase generado por una sonda de amplitud es una componente cruzada pequena. Su precision relativa tardia no esta acreditada; no debe presentarse como un resultado cuantitativo resuelto. La amplitud, la corriente y la fuerza principales si concuerdan entre los dos pasos de tiempo. Esto permite continuar el desarrollo pertinente al sistema final sin repetir esta bateria solo para resolver una cola pequena. El estado historico sigue siendo TEMPORAL_REFINEMENT_NOT_MET.

## Magnitudes y figuras

- physical_case: Fixed bath T=0.9 K, Tc=8.65 K; 65x65 graph vortex-core control on square [-6,6] ell0 in x and y, spacing 0.1875 ell0 and 256 fixed boundary nodes; 256 Matsubara terms; no photon, no strip experiment and no circuit trajectory.
- baseline: d_b(t) is the unperturbed evolving core. It is not stationary and is subtracted at the same time.
- amplitude_probe: Initial d_A-d_b = 0.001 d0 b(r), b=max(0,1-r^2/R^2)^3; R=4 ell0 from admitted operator plan.
- angular_phase_probe: Initial d_P-d_b = 0.001 i d0 (x/R) b(r). Cartesian perturbation, not a global phase shift.
- node_norm: ||z||_M=sqrt(sum_i m_i |z_i|^2), including fixed nodes, m_i dimensionless cell area.
- current_norm: ||j||=sqrt(sum_edges_positive_c j_e^2/c_e), c_e graph conductance in dimensionless action.
- displacement: d=Delta/(k_B Tc). Curve: ||d_probe(t)-d_b(t)||_M / ||d_probe(0)-d_b(0)||_M; unitless, not a point amplitude.
- force_density: q_i=G_i/m_i, where G is discrete thermal free-energy gradient. Curve uses probe-minus-base force; dimensionless variational density, not newtons.
- phase_torque_density: tau_i=Im(conj(d_i) G_i)/m_i. Curve uses tau_probe-tau_base; dimensionless force for condensate phase, not mechanical torque.
- current: Derivative of the same dimensionless spectral action with respect to edge gauge link. Probe-minus-base graph current; no ampere conversion is asserted.
- time: Real physical time in ps; tau=t/[hbar/(2 k_B Tc)] is internal solver time only.
- refinement: Primary minus refined trajectory in the same observable and norm. Original gate: error <= 1e-8 + 0.005 max(initial norm, remaining refined norm).
- contextual_scale: Maximum initial amplitude/phase signal for the SAME observable. A diagnostic contextualization, not a replacement gate.
- free_energy: Change in dimensionless finite-Matsubara thermal free-energy action from d0 at fixed bath. It is not conserved internal energy.
- plots_must_distinguish: Response magnitude, refinement difference, and remaining-signal denominator. State which field, probe, norm, normalization, units and baseline appear in each caption.

Verificados 71 archivos NPZ y 55 pasos aceptados. Reproduccion independiente desde campos: diferencia maxima 2.39e-18. No se ejecutaron nuevas ecuaciones fisicas.
