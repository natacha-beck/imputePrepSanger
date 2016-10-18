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
    && mkdir ressources \
    && mkdir ressources/data \
    && mkdir ressources/HRC_refSites \
    && mkdir ressources/strand

#RUN cd /imputePrepSanger/ressources/data/ && echo | ls

# These files are from the git clone, and need to move into specific directories
RUN mv HRC-1000G-check-bim_modified.pl ressources/HRC_refSites/
RUN mv ucsc2ensembl.txt ressources/HRC_refSites/
RUN mv update_build.sh ressources/strand/

WORKDIR tools
RUN wget  https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 \
    && bunzip2 bcftools-1.3.1.tar.bz2 \
    && tar -xvf bcftools-1.3.1.tar \
    && cd bcftools-1.3.1 \
    && mkdir bin \
    && make \
    && make install


RUN mv bcftools-1.3.1/bcftools bcftools-1.3.1/bin

#Now we install plink 1.9
WORKDIR plink/
RUN wget https://www.cog-genomics.org/static/bin/plink161010/plink_linux_x86_64.zip \
    &&  unzip -a plink_linux_x86_64.zip

# Now we need some more data
#WORKDIR ../../ressources/HRC_refSites/
#RUN wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
#    && gunzip -f human_g1k_v37.fasta.gz

#RUN wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1/HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz
#    && gunzip -f HRC.r1.GRCh37.autosomes.mac5.sites.tab.gz

WORKDIR ../../
