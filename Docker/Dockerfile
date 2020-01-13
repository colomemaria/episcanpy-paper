
FROM jupyter/minimal-notebook

USER root

USER $NB_UID

# Install Python 3 packages
RUN conda install --quiet --yes \
      --channel bioconda \
      --channel conda-forge \
    'matplotlib-base=3.1.*' \
    'numba=0.45*' \
    'pandas=0.24*' \
    'scikit-learn=0.21*' \
    'scipy=1.3*' \
    'seaborn=0.9*' \
    'statsmodels=0.10*' \
    'python-igraph==0.7*' \
    'louvain==0.6*' \
    'pysam' \
    && \
    conda clean --all -f -y && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER


RUN pip install --upgrade --no-cache-dir --upgrade-strategy=only-if-needed \
    cmake==3.13.3 && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

RUN pip install --upgrade --no-cache-dir --upgrade-strategy=only-if-needed \
    MulticoreTSNE==0.1 \
    numpy==1.16.4 \
    scanpy==1.4.4.post1 && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

RUN pip install git+https://github.com/colomemaria/epiScanpy@60484de8ccb44cf53ead1cafece11cf5e8e9dd42 && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

