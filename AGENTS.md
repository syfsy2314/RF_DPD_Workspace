# Repository Guidelines

## Project Structure & Module Organization

This is a MATLAB workspace for RF power-amplifier modeling and digital predistortion (DPD) experiments.

- `common/` contains reusable code: PA models under `models/`, signal and metric helpers under `tools/`, and shared experiment utilities under `experiments/`.
- `projects/` contains algorithm-specific work, including `LD-PLS`, `MG`, `MR`, `SBL`, and `MP_DPD`. Keep project entry scripts and project-only helpers together.
- `datasets/` is reserved for shared input data. Small project-specific `.mat` files may live in a project's `data/` directory.
- `.matlab/startup.m` implements path initialization; root `startup.m` is the MATLAB-discoverable entry point and must remain lightweight.

Open the repository root in VS Code so MATLAB starts in the correct folder and recursively adds the workspace to its path.

## Build, Test, and Development Commands

There is no compilation step. Run commands from the repository root:

```powershell
matlab -batch "startup; which get_NMSE"
matlab -batch "checkcode('common/tools/metrics/get_NMSE.m')"
matlab -batch "results=runtests('tests'); assertSuccess(results)"
```

The first command verifies startup and path resolution. `checkcode` performs MATLAB static analysis. The final command runs the test suite once `tests/` exists. Run experiments with `run('projects/<name>/main.m')`; confirm required data and toolboxes first.

## Coding Style & Naming Conventions

Use four-space indentation, one statement per line, and spaces around binary operators. A function's name must match its filename. Follow existing descriptive prefixes such as `get_`, `plot_`, and `create_`; use `main.m` only for project entry scripts. Add `%` comments for algorithm assumptions and units, not line-by-line narration. Avoid new hard-coded machine paths. Resolve files with `fullfile` relative to the project or script location. 不用过多防御性编程，代码以简洁为主；仅在存在明确风险或公共接口确实需要时增加输入校验和异常分支。

## Testing Guidelines

Add MATLAB unit tests under each project's `tests/` directory, such as `projects/Alternating Iterations/tests/`. Name files `test_<subject>.m` and test methods by behavior, such as `test_identicalSignalsHaveZeroNMSE`. Use deterministic inputs and fixed random seeds. Cover nominal behavior, dimensions, complex-valued signals, and invalid arguments. Keep generated results in the same project's `results/` directory, and do not overwrite committed `.mat` fixtures during tests.

## Commit & Pull Request Guidelines

No commit history is currently available, so use concise imperative messages, preferably scoped: `feat(ld-pls): add regularization sweep` or `fix(metrics): handle zero reference power`. Keep commits focused. Pull requests should describe the algorithmic change, MATLAB release and required toolboxes, data provenance, validation commands, and numerical before/after results. Include plots when signal behavior changes, but exclude generated figures and large datasets unless explicitly required.

## Agent Connection Handling

When a network or remote-tool connection fails, reconnect and retry automatically up to three times. If all three attempts fail, stop retrying and report exactly `网络连接问题` to the user.
