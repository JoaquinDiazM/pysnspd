# Repository working policy

- Work directly on `main` for routine `pysnspd` tasks.
- Do not create a feature branch or pull request unless the user explicitly asks for one.

## Geminga computation policy

- Use the non-administrator `jdiaz` account. No `sudo` or administrator access.
- Do not launch computations known or reasonably expected to take more than five minutes.
- Write those commands, their purpose, expected outputs and estimated resources to
  `/home/jdiaz/GEMINGA_COMMANDS.md`. Preserve the existing entries.
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
