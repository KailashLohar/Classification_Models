# Base image
FROM --platform=linux/arm64/v8 ubuntu:22.04

# Set Label
LABEL maintainer="Kailash Lohar <kailashlohariitkgp@gmail.com>"

# Set Working Directory
WORKDIR /home

# Install necessary utilities
RUN apt-get update && apt-get install -y wget curl git && \
    wget -O miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh && \
    bash miniconda.sh -b -p /miniconda && rm miniconda.sh && rm -rf /var/lib/apt/lists/*

# Set Conda in the environment path
ENV PATH /miniconda/bin:$PATH

# Initialize conda
RUN conda init bash

# Copy your FastAPI YAML environment file into the container
COPY fastapi.yaml .

# Create the environment using your YAML file
RUN conda env create -f fastapi.yaml && conda clean --all --yes

# Expose FastAPI port (8502)
EXPOSE 8502

# Copy SSL certificates into the container
COPY cert.pem /home/cert.pem
COPY key.pem /home/key.pem

# Start FastAPI via Uvicorn using SSL
CMD ["/bin/bash", "-c", "source /miniconda/etc/profile.d/conda.sh && conda activate fastapi_env && uvicorn app:app --host 0.0.0.0 --port 8502 --ssl-keyfile /home/key.pem --ssl-certfile /home/cert.pem"]
