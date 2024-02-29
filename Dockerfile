FROM continuumio/miniconda3:4.7.12
LABEL authors="Bayram Boukhari"

RUN apt-get update --fix-missing




ENV PATH /Bracken:$PATH

COPY environment.yaml .
RUN conda env update -n root -f environment.yaml && \
    conda clean -afy && pip install git+https://github.com/exels/MultiQC.git