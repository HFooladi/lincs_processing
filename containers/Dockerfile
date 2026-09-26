# One image for every Nextflow process (target `pipeline`) and a slightly
# larger one for the FastAPI/Gradio app (target `serve`).
#
#   docker build -f containers/Dockerfile --target pipeline -t lincs-processing:0.3.0 .
#   docker build -f containers/Dockerfile --target serve    -t lincs-processing-serve:0.3.0 .
#
# PyTorch wheels from PyPI bundle their CUDA libraries, so the base image only
# needs a working NVIDIA driver on the host (docker --gpus all).

FROM python:3.13-slim-bookworm AS pipeline

# procps: Nextflow uses `ps` to collect task metrics.
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps libxrender1 libxext6 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:0.9 /uv /uvx /bin/

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    PATH=/opt/venv/bin:$PATH \
    MPLCONFIGDIR=/tmp

WORKDIR /app

# Dependencies first so code changes do not rebuild the torch layer.
COPY pyproject.toml uv.lock README.md ./
RUN uv sync --frozen --no-dev --group pipeline --no-install-project

COPY src ./src
RUN uv sync --frozen --no-dev --group pipeline --no-editable \
    && lincs-pipeline --help > /dev/null

FROM pipeline AS serve

RUN uv sync --frozen --no-dev --group pipeline --group serve --no-editable
COPY serve ./serve
EXPOSE 8000
CMD ["uvicorn", "serve.app:create_app", "--factory", "--host", "0.0.0.0", "--port", "8000"]
