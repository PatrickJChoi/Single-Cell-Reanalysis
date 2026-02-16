FROM condaforge/mambaforge:latest

WORKDIR /app

# Create conda environment from environment.yml
COPY environment.yml .
RUN mamba env create -f environment.yml && mamba clean --all -y

# Activate environment by default
ENV PATH=/opt/conda/envs/scr_smoke/bin:$PATH
ENV CONDA_DEFAULT_ENV=scr_smoke

# Copy project files
COPY . .

# Create output directories
RUN mkdir -p data/raw outputs reports results/tables results/figures

# Generate toy data so smoke test is ready to run
RUN python src/pipeline/make_toy_data.py

# Default: run the smoke test
CMD ["python", "src/pipeline/run.py", "--config", "configs/smoke.yaml"]
