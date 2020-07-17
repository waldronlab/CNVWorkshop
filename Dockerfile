FROM bioconductor/bioconductor_docker:bioc2020

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); BiocManager::install('waldronlab/CNVWorkshop', update = TRUE, ask = FALSE)"
