# Contributing

Thank you for helping improve this workflow. Bug reports, questions, documentation fixes and code contributions are all welcome.

## Reporting a problem

Open an [issue](../../issues/new/choose) and include:

- the rule that failed and its log from `results/logs/`
- the output of `snakemake --version` and `conda list -n yfacet-snakemake | head`
- your operating system and whether you run on SLURM or locally

## Suggesting a change

1. Fork the repository, or create a branch if you are a collaborator: `git checkout -b fix/short-description`.
2. Make the change and keep it focused on one thing.
3. Run the checks:
   ```bash
   python tests/check_workflow.py
   bash tests/run_toy_tests.sh
   ```
4. Update the documentation in `docs/` if behaviour or parameters change.
5. Add a line to `CHANGELOG.md` under "Unreleased".
6. Open a pull request and describe what changed and why.

## Code style

- **One step, one script.** Each script in `workflow/scripts/` does one job, has a header comment explaining it, and takes its inputs as command-line options (`-h` must work).
- **Comment with `#`.** Explain *why* a step is done, not only what it does.
- **No hard-coded values.** Thresholds and parameters belong in `config/config.yaml`.
- **Report numbers.** A step that produces a number the manuscript needs writes it to a `*.metrics.tsv` file (key and value), so it reaches `results/summary_numbers.tsv`.
- **Environments.** New tools go into the matching file in `workflow/envs/`, or a new one. Also add the environment name to `ENV_NAMES` in `workflow/rules/common.smk`.

## Scientific changes

Changes that alter results, such as thresholds, filters or methods, should be discussed in an issue first. They must be reported in the manuscript Methods and listed in `docs/methods-changes.md`.

## Code of Conduct

By participating you agree to follow the [Code of Conduct](CODE_OF_CONDUCT.md).
