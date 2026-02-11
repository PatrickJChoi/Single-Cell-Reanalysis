\# Contributing



Thanks for your interest in contributing!



\## Quick guidelines



\- Open an issue for bugs, questions, or enhancement ideas.

\- Keep changes focused and small when possible.

\- For analysis changes, prefer config-driven behavior (see `configs/`).



\## Development setup



1\) Create the environment:



&nbsp;   conda env create -f environment.yml



2\) Activate it:



&nbsp;   conda activate scr\_smoke



\## Before you open a PR



Run the smoke test locally:



&nbsp;   python src/pipeline/make\_toy\_data.py

&nbsp;   python src/pipeline/run.py --config configs/smoke.yaml



Make sure it completes successfully and produces outputs in `reports/smoke/`.



\## Notes



\- Generated data/figures are not committed (see `.gitignore`).

\- Real datasets are not included; see `docs/Data\_Acquisition.md`.



