FROM quay.io/qiime2/amplicon:2026.1

ENV DEBIAN_FRONTEND=noninteractive

# System deps
RUN apt-get update && apt-get install -y \
    curl git build-essential sudo \
    && rm -rf /var/lib/apt/lists/*

# Node (for mystmd)
RUN curl -fsSL https://deb.nodesource.com/setup_20.x | bash - \
    && apt-get install -y nodejs

RUN npm install -g mystmd@1.8.2

# Python dev tools inside base env
RUN /opt/conda/bin/pip install ipython pdbpp rich

# vscode user
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=1000

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# 👇 KEY FIX: make conda usable in interactive shells
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> /home/$USERNAME/.bashrc

# Install debug tools and dependencies
RUN apt-get update && apt-get install -y gdb cmake zlib1g-dev

USER $USERNAME
WORKDIR /workspace

RUN conda run -n qiime2-amplicon-2026.1 Rscript -e "install.packages('languageserver', repos='https://cloud.r-project.org')"


# IMPORTANT: use login bash
ENTRYPOINT ["bash", "-l"]