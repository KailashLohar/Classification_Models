FROM --platform=linux/arm64/v8 ubuntu:22.04

LABEL maintainer="Kailash Lohar <kailashlohariitkgp@gmail.com>"

WORKDIR /home

# Install basic utilities and Miniconda
RUN apt-get update && apt-get install -y wget curl git python3 python3-pip && \
    wget -O miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh && \
    bash miniconda.sh -b -p /miniconda && rm miniconda.sh && rm -rf /var/lib/apt/lists/*

# Set PATH for conda
ENV PATH /miniconda/bin:$PATH

# Configure conda channels for optimal dependency resolution
RUN echo "channels:" > /miniconda/.condarc && \
    echo "  - conda-forge" >> /miniconda/.condarc && \
    echo "channel_priority: strict" >> /miniconda/.condarc

# Copy environment configuration file
COPY streamlit.yaml .

# Remove existing environments and install dependencies
RUN conda env remove -n streamlit_env || true && \
    conda env create -f streamlit.yaml || (sleep 5 && conda env create -f streamlit.yaml) && \
    conda clean --all --yes

# Expose Streamlit port
EXPOSE 8501

# Default command to run Streamlit app
CMD ["/bin/bash", "-c", "source /miniconda/etc/profile.d/conda.sh && conda activate streamlit_env && streamlit run main.py --server.port 8501"]
