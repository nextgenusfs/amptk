# start with miniconda3 as build environment
FROM continuumio/miniconda3 AS build

# Update, install mamba and conda-pack:
RUN conda update -n base -c defaults --yes conda && \
    conda install -c conda-forge -n base --yes mamba conda-pack

# Install amptk deps from bioconda
# here specifying specific versions to be able to set ENV below
RUN mamba create -c conda-forge -c bioconda -c defaults \
    -n amptk --yes "python>=3.6" biopython distro psutil \
    "python-edlib>=1.2.1" "vsearch>=2.15.0" natsort pandas \
    matplotlib-base seaborn biom-format r-base "bioconductor-dada2>=1.12.1" \
    bioconductor-phyloseq r-plotly r-htmltools r-dt pigz pyfastx \
    mafft fasttree requests minimap2 && conda clean -a -y

# Since we want the most recent, install from repo
SHELL ["conda", "run", "-n", "amptk", "/bin/bash", "-c"]
RUN python -m pip install git+https://github.com/nextgenusfs/amptk.git

# package with conda-pack
RUN conda-pack -n amptk -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

# We've put venv in same path it'll be in final image
RUN /venv/bin/conda-unpack

# Now build environment
FROM debian:buster AS runtime

# Copy /venv from the previous stage:
COPY --from=build /venv /venv

# add it to the PATH and add env variables
ENV PATH="/venv/bin:$PATH"

# When image is run, run the code with the environment
SHELL ["/bin/bash", "-c"]
CMD amptk
