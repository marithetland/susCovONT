#! /bin/sh

#A wrapper script to install all of the necessary tools and conda env has been made in scripts/install.sh. The only thing you need to set it INSTALL_DIR which is where the repositories will be installed.
#INSTALL_DIR="/media/marit/Programs/"  #Change to your install di
#cd $INSTALL_DIR
#git clone https://github.com/marithetland/susCovONT  #Clone this repo
#cd susCovONT
#bash /install.sh $INSTALL_DIR #Install all necessary tools
#This will install the necessary repos and conda environment and will update the config file scripts/config.cfg to have the correct paths. The only other thing you should do is change the number of cpus to use in

## Usage
#bash susCovONT/scripts/install.sh /path/to/install_dir

echo "##### Checking input ..."
## Get installation directory from user input and check that the directory exists, then enter dir
if [ -d $1 ]; then
    echo "Found $1."
    INSTALL_DIR=$(echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")")
    cd $INSTALL_DIR
else
    echo "Error: Directory $1 does not exist. Please provide a path that exists."
    exit
fi

# Check that this repo exists under INSTALL_DIR
if [ -d ${INSTALL_DIR}/susCovONT ]; then
    echo "Found ${INSTALL_DIR}/susCovONT."
else
    echo "Error: Directory ${INSTALL_DIR}/susCovONT does not exist. Please make sure you have cloned the susCovONT repo in the INSTALL_DIR folder."
    echo "cd ${INSTALL_DIR} ; git clone https://github.com/marithetland/susCovONT"
    exit
fi

# Check that conda and docker is installed
testcmd () {
    command -v "$1" >/dev/null
}

if testcmd conda ; then
    echo "Found conda."
else
    echo "Did not find conda in path. Have you installed it?"
    exit
fi

if testcmd docker ; then
    echo "Found docker."
else
    echo "Did not find docker in path. Have you installed it?"
    exit
fi

## Check with user
echo "##### The following will be installed: " 
echo "- Artic nextflow pipeline in ${INSTALL_DIR}/ncov2019-artic-nf"
echo "- Conda environment for artic nextflow pipeline in ${INSTALL_DIR}/conda_for_covid/"
echo "- Pangolin (and pangolin environment) for artic nextflow pipeline in ${INSTALL_DIR}/pangolin"
echo "- Nextstrain/nextclade image will be pulled"
echo "- Config file for susCovONT will be updated with the paths above."

Read -p "#### Do you wish to proceed the installation with these (y/n)?" choice
case "$choice" in 
  y|Y ) echo "Excellent - installing now ...";;
  n|N ) echo "Aborting - you did say no ..."; exit;;
  * ) echo "Aborting - invalid character. Please respond to prompt with 'y' or 'Y'."; exit;;
esac


## Git clone artic nextflow pipeline
echo "##### Cloning artic nextflow pipeline to ${INSTALL_DIR}/ncov2019-artic-nf"
cd $INSTALL_DIR
git clone https://github.com/connor-lab/ncov2019-artic-nf
sed -i.bak 's/artic=1.1.3/artic=1.2.1/g' ${INSTALL_DIR}/ncov2019-artic-nf/environments/nanopore/environment.yml #Change "1.2.1" to desired artic version
cp ${INSTALL_DIR}/susCovONT/scripts/articQC.py ncov2019-artic-nf/bin/qc.py #Update the QC script with the one from this repository.

## Create conda env to use as --cache for nextflow pipeline
echo "##### Creating conda environment with artic v1.2.1 (uses conda)"
echo "##### Bear with us, this takes a little while..."
cd $INSTALL_DIR
CONDA_LOCATION=$(echo "${INSTALL_DIR}/conda_for_covid/") #Change this to match your system
mkdir $CONDA_LOCATION
conda env create --prefix ${CONDA_LOCATION}/work/conda/artic-2c6f8ebeb615d37ee3372e543ec21891 --file ${INSTALL_DIR}/ncov2019-artic-nf/environments/nanopore/environment.yml

## Download primer schemes to be used in offline mode
echo "##### Downloading primer schemes for offline mode to ${INSTALL_DIR}/primer-schemes"
cd $INSTALL_DIR
git clone https://github.com/artic-network/primer-schemes.git
##Check that folder didn't exist and that it now does

## Install pangolin
echo "##### Installing pangolin (uses conda)"
cd $INSTALL_DIR
git clone https://github.com/cov-lineages/pangolin.git
cd ${INSTALL_DIR}/pangolin
conda env create -f environment.yml
bash -c "source activate pangolin ; conda activate pangolin ; python setup.py install ;conda deactivate"

## Pull nextclade image
echo "##### Pulling nextclade (uses docker)"
cd $INSTALL_DIR
docker pull nextstrain/nextclade

## Change config file
echo "##### Updating paths in file ${INSTALL_DIR}/susCovONT/scripts/config.cfg"
cd $INSTALL_DIR
sed -i "s|nf_location = /home/marit/Programs/ncov2019-artic-nf/|nf_location = ${INSTALL_DIR}/ncov2019-artic-nf/|" ${INSTALL_DIR}/susCovONT/scripts/config.cfg
sed -i "s|conda_location = /home/marit/Programs/conda_for_covid/work/conda/|conda_location = ${INSTALL_DIR}/conda_for_covid/work/conda/|" ${INSTALL_DIR}/susCovONT/scripts/config.cfg
sed -i "s|schemeRepoURL = /home/marit/Programs/primer-schemes/|schemeRepoURL = ${INSTALL_DIR}/primer-schemes/|" ${INSTALL_DIR}/susCovONT/scripts/config.cfg

#Done
