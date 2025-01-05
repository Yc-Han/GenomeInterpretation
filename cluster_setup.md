# Installation Guide on Cluster

Copyright \@ Xiao-Yin To

## 1. Install Miniconda and TensorFlow

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod 755 Miniconda3-latest-Linux-x86_64.sh 
bash ./Miniconda3-latest-Linux-x86_64.sh
# press enter until it asks you to say 'yes'
# afterwards ENTER for default location, no to second question

echo '. ~/miniconda3/etc/profile.d/conda.sh' >> ~/.bashrc
source ~/.bashrc # or restart shell

conda create --name tf python=3.9 # must-be
conda activate tf
conda install -c conda-forge cudatoolkit=11.8.0
pip install nvidia-cudnn-cu11==8.6.0.163
pip install --upgrade pip

pip install tensorflow==2.13.* # must-be
  
conda install r-base=4.4.2
conda install r-devtools
conda install h5py
```

## 2. Install DeepG

```R
install.packages("tensorflow")
tensorflow::install_tensorflow()
devtools::install_github("GenomeNet/deepG")
```

In some cases, `tensorflow::install_tensorflow()` does not work.

Then, `reticulate::py_install("tensorflow=2.13.*", pip=TRUE)` can be used as an alternative.

## 3. Configure the system paths

```bash
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib/' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
```
