FROM ubuntu:16.04
LABEL maintainer="frouin.eleonore@gmail.com"

#### Install Essential
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update \
    && apt-get -y install pigz \
                          libcurl4-gnutls-dev \
                          openjdk-8-jdk \
                          libssl-dev \
                          libxml2-dev


#### Install MiniConda
RUN apt-get -qq update && apt-get -qq -y install curl bzip2 \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && apt-get -qq -y remove curl bzip2 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes
ENV PATH /opt/conda/bin:$PATH


#### Install Pandaseq
RUN conda install -c bioconda pandaseq=2.11


#### Install Trimmomatic
RUN conda install -c bioconda trimmomatic=0.36


#### Install bwa
RUN conda install -c bioconda bwa=0.7.12


#### Install samtools
RUN conda install -c bioconda samtools=1.3.1


#### Install prokka and prodigal
ENV PERL5LIB=/usr/local/lib/perl5/site_perl/5.22.0/
RUN apt-get -qq update \
	&& apt-get -qq -y install less 
RUN conda install -c conda-forge -c bioconda prokka=1.11


#### Install Diamond
RUN conda install -c bioconda diamond=0.8.36


#### Install Blast and RPSblast
RUN conda install  -c bioconda blast=2.2.31

#### Install bedtools
RUN conda install -c bioconda bedtools=2.26.0


#### Install make and cie.
RUN apt-get update && apt-get -y install build-essential


#### Install cd-hit
ADD https://github.com/weizhongli/cdhit/releases/download/V4.6.6/cd-hit-v4.6.6-2016-0711.tar.gz /tmp/cd-hit.tar.gz
RUN mkdir /tmp/cd-hit \
    && tar xzf /tmp/cd-hit.tar.gz --directory /tmp/cd-hit --strip-components=1
RUN cd /tmp/cd-hit/cd-hit-auxtools \
    && make 
RUN mv /tmp/cd-hit/cd-hit-auxtools/ /opt/cd-hit-auxtools/ \
    &&  ln -s /opt/cd-hit-auxtools/cd-hit-dup /usr/bin/cd-hit-dup

#### Install IDBA-UD
ADD https://github.com/loneknightpy/idba/releases/download/1.1.2/idba-1.1.2.tar.gz /tmp/idba.tar.gz
RUN mkdir /tmp/idba
RUN tar xzf /tmp/idba.tar.gz --directory /tmp/idba --strip-components=1 \
    && sed --in-place 's/kMaxShortSequence = 128;/kMaxShortSequence = 1024;/' /tmp/idba/src/sequence/short_sequence.h
RUN cd /tmp/idba \
	&& ./configure  \
	&& make  \
	&& make install 
RUN mv /tmp/idba/bin/* /usr/local/bin/


#### Install snakemake
RUN conda install -c bioconda -c conda-forge snakemake=3.5.5


#### Install python2.7 and modules
RUN apt-get -qq update \
	&& apt-get -qq -y install python-pip \
	&& python2.7 -m pip install biopython \
	&& python2.7 -m pip install pandas 


#### Clean tmp dir
RUN rm -fr /tmp/*

#### Create TreeView
RUN mkdir /root/Data 
ADD Snakefile /root
ADD config.yaml /root
ADD databases /root/databases
ADD scripts /root/scripts
RUN ln -s /root/scripts/* /usr/bin/







