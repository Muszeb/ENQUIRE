: '
A) CONFIGURE INSTALLATION
Change the following variable to customize the installation path:
'

mambapath=$HOME/bin

: '
B) INSTALL
To install ENQUIRE using a virtual environment, we need to 
  1) Install a virtual environment manager (here, we use micromamba - https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html);
  2) Create and activate virtual environment called ENQUIRE;  
  3) Install Python, R and all the required library dependencies;
  4) Install EDirect - https://www.ncbi.nlm.nih.gov/books/NBK179288/ 
'

# minimal dependencies 
sudo apt-get update
sudo apt-get install -y curl

# Download micromamba and export path to micromamba (possibly update .bashrc)
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj "${mambapath}/micromamba"
export PATH="${mambapath}:$PATH"
echo "export PATH=${mambapath}:$PATH" >> ${HOME}/.bashrc # not necessary but recommended

# This command allows to non-interactively run micromamba commands
#eval "micromamba shell hook -s posix --shell bash"
eval "$(micromamba shell hook --shell bash)"

# Create an ad hoc environment.
# If you have another virtual environment manager installed (e.g miniconda), the environment might be installed under '$HOME/othermanager/envs'
micromamba create -n ENQUIRE 
micromamba activate ENQUIRE

# Install EDirect under your HOME directory - manually adding the edirect path to .bashrc keeps the latter cleaner
yes n | sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
echo "export PATH=${HOME}/edirect:${PATH}" >> "${HOME}/.bashrc" # necessary 

# Install the python packages and R under the ENQUIRE environment, from the ENQUIRE.yml file
micromamba install -y -q -f input/ENQUIRE.yml
pip install -r input/PackagesNotFound.txt # these packages can't be installed by environment managers

# Install R libraries under the ENQUIRE environment (remember to "activate ENQUIRE"!)

$(which Rscript) code/install_R_libraries.R 

# Install Pandoc 

sudo apt-get install -y pandoc

# Clean environment

micromamba clean --all --yes