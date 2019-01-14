FROM ubuntu:18.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV MPLBACKEND TkAgg

RUN apt-get update --fix-missing && \
    apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    curl \
    build-essential \
    git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]

RUN conda update conda
RUN conda config --add channels defaults
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install --quiet --yes scanpy pip dask dca pypairs numba multicore-tsne cython joblib pandas numpy scipy matplotlib
RUN conda clean -tipsy

RUN pip --no-cache-dir install phate fa2 bbknn magic-impute leidenalg vpolo annoy scprep zappy bioservices

#RUN git clone --recursive git://github.com/KrishnaswamyLab/PHATE.git
#RUN cd PHATE/Python
#RUN python setup.py install

RUN pip install git+https://github.com/flatberg/scanpy-scripts.git