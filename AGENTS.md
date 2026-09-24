# Repository working policy

- Work directly on `main` for routine `pysnspd` tasks.
- Do not create a feature branch or pull request unless the user explicitly asks for one.

## Geminga computation policy

- Use the non-administrator `jdiaz` account. No `sudo` or administrator access.
- Parallelize independent cases and expensive independent spectral queries.
  Recheck CPU topology, effective affinity/quotas and available memory at runtime.
  Allocate at most 90% of the available computational resources; reserve complete
  physical cores where possible. The 2026-09-23 inventory is 16 physical cores,
  32 logical CPUs and about 123 GiB RAM: at most 28 logical CPUs, leaving two
  complete cores free. This is a ceiling, not a requirement to fill idle slots.
- Count nested workers and numerical-library threads within one shared budget.
  Pin the task to the admitted CPU set, limit BLAS/OpenMP threads per worker,
  and reduce concurrency when memory or a quota is tighter. Report the actual
  worker count and progress/ETA. Never assume that independent cases need to run
  successively or that an untested parallel speedup makes a long job lightweight.
- Do not launch computations known or reasonably expected to take more than five minutes.
- Write those commands, their purpose, expected outputs and estimated resources to
  `/home/jdiaz/GEMINGA_COMMANDS.md`. Preserve the existing entries.
- Also include the exact long-run command explicitly in the chat, in a copyable
  code block. Do not refer the user only to the command notebook.
- Leave useful commands the user may want to try in that same file, including
  reproducible lightweight diagnostics.
- For a required long calculation, stop the dependent work and wait for the user
  to run it and return the terminal output with `reinicia` or an equivalent request.
  Do not use polling loops, scheduled wakeups or background jobs to bypass this handoff.
- Use a bounded timeout for lightweight checks with uncertain runtime. A timed-out
  check is incomplete; record the unrestricted user-run command instead of restarting
  it repeatedly or splitting a long calculation to evade the five-minute limit.
- The tag `v1.0.0` freezes the thesis implementation with model 0.4 documentation.
  Later experimental implementation belongs in subsequent commits, not in that tag.

## Practical numerical acceptance

- Prefer revising an unnecessarily strict tolerance to postponing a small
  discrepancy. Tie the revised margin to the observable, its signal scale and
  the purpose of the calculation; record both the original result and the new
  acceptance decision. A tiny residual signal need not block development.
- Do not use tolerance changes to conceal a large physical error, nonfinite
  states, instability, missing coupling terms or material conservation defects.
  Explain such failures and repair the cause before promoting that capability.
- Reuse published and inherited algorithms before adding a new solver. Prefer
  the thesis/pyTDGL KWT local Euler update for the next temporal comparison.
  Its quadratic amplitude solve does not make it second-order in time; establish
  accuracy by step refinement of the relevant observables.
- Reports must label each plotted quantity, units, reference subtraction, norm,
  normalization and physical case explicitly.

## Initial photon-transient observation window

- The first target is relative detection latency for the 80 nm Korzh (2020)
  device. Evolve through the agreed `V_out` trigger plus a declared extra window;
  full nanosecond-scale recovery is not the initial target.
- This changes only the observation/integration horizon. Preserve the same
  physical equations, material and device parameters, deposition, mesh, boundary
  conditions, full thesis circuit and numerical accuracy policy when comparing
  short and extended windows. Do not shorten circuit time constants or discard
  slow physical states to make the window shorter.
- Record a first crossing from accepted time steps with an identified voltage
  observable, threshold and crossing confirmation; record the extra margin.
  Do not invent these physical measurement choices. An event that never triggers
  before its finite time ceiling is censored, not assigned a fictitious latency.
- The existing production `latency` early-stop mode waits for a peak; it is not
  yet the requested trigger-plus-margin semantics. Production changes require
  the later implementation step, not a silent change in these static controls.
