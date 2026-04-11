FROM quay.io/qiime2/amplicon:2026.1

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    git \
    build-essential \
    sudo \
    && rm -rf /var/lib/apt/lists/*

# ---- Create vscode user ----
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# ---- Install Node.js (as root) ----
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - \
    && apt-get install -y nodejs

# ---- Install mystmd globally ----
RUN npm install -g mystmd

# ---- Python dev/debug tools (inside QIIME env later) ----
RUN /opt/conda/bin/pip install \
    ipython \
    pdbpp \
    rich

# ---- Switch to vscode user ----
USER $USERNAME
WORKDIR /workspace

# Activate QIIME2 automatically
RUN echo "conda activate qiime2-amplicon-2026.1" >> ~/.bashrc

ENTRYPOINT "/bin/bash"