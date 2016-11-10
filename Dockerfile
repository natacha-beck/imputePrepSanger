FROM centos

MAINTAINER Marie Forest <marie.forest@ladydavis.ca>

# Install prerequisite
RUN yum update -y

# Install basic packages
RUN yum install -y bzip2 \
                   gcc \
                   git \
                   make \
                   perl \
                   unzip \
                   wget \
                   zlib-devel

RUN git clone https://github.com/eauforest/imputePrepSanger.git

WORKDIR /imputePrepSanger/
#RUN  echo | ls

RUN mkdir results \
    && mkdir tools \
    && mkdir tools/plink \
    && mkdir ressources \
    && mkdir ressources/data \
    && mkdir ressources/HRC_refSites \
    && mkdir ressources/strand

#RUN cd /imputePrepSanger/ressources/data/ && echo | ls

# These files are from the git clone, and need to move into specific directories
RUN mv HRC-1000G-check-bim_modified.pl ressources/HRC_refSites/
RUN mv ucsc2ensembl.txt ressources/HRC_refSites/
RUN mv update_build.sh ressources/strand/
RUN mv bcftools-1.3.1.tar.bz2 tools/
RUN mv plink_linux_x86_64.zip tools/plink/

WORKDIR tools
RUN  bunzip2 bcftools-1.3.1.tar.bz2 \
    && tar -xvf bcftools-1.3.1.tar \
    && cd bcftools-1.3.1 \
    && mkdir bin \
    && make \
    && make install


RUN mv bcftools-1.3.1/bcftools bcftools-1.3.1/bin

#Now we unzip plink 1.9
WORKDIR plink/
RUN  unzip -a plink_linux_x86_64.zip


WORKDIR ../../
