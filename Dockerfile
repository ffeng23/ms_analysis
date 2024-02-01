FROM rocker/rstudio:4.2.2
#create a folder in container and entering it.
#something it is not working well at the root directory.
WORKDIR /home/rstudio/clinicalTrial
#note about directory, we use the bind volume to the container volum main

RUN sudo apt-get update && apt-get install -y --no-install-recommends libxt6 \
    && apt-get install -y libz-dev sudo cmake \
    && apt-get install -y libpng-dev libfontconfig1-dev \
    && apt-get install -y libharfbuzz-dev libfribidi-dev \
    && apt-get install -y libfreetype6-dev libtiff5-dev libjpeg-dev \
    && apt-get install -y libxml2-dev 

RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# approach two
# to set up the 
COPY renv.lock renv.lock
#COPY renv.lock /home/renv.lock
RUN  mkdir -p renv
RUN  mkdir -p R_code
RUN  mkdir -p output_flow
RUN  mkdir -p Data
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
#COPY renv/settings.json renv/settings.json
COPY renv/settings.dcf renv/settings.dcf
COPY R_code/* R_code/

# Note for above, we don't have renv/settings.json somehow (?),
# it is not necessary!!! we could just need .Rprofile and activate.R and renv.lock
#

#change owner ship to rstudio, RECURSIVELY with -R
# and then restor the renv library as the user "rstudio"
RUN chown -R rstudio . \
    && sudo -u rstudio R -e "renv::restore()"
#RUN chown rstudio /home/renv.lock

RUN sudo -u rstudio echo "setwd(\"/home/rstudio/clinicalTrial/\")" > /home/rstudio/.Rprofile

RUN sudo -u rstudio echo "renv::load()" >> /home/rstudio/.Rprofile 
#    && echo ".rs.restartR()" >>/home/rstudio/.Rprofile 


