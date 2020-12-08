FROM nfcore/base:1.9
LABEL authors="Chela James" \
      description="Docker image containing all software requirements for the umi_preprocessing"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-sarek-2.6.1/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name umi-umi_preprocessing-0.0.1 > umi-umi_preprocessing-0.0.1.yml