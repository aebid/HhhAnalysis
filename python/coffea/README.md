# Instruction to setup and run coffea framework

## Set up conda virtual environment on lxplus
Coffea framework usually is executed under virtual environment that would not intefere with other settings. To set up 
the virtual environment on lxplus
  - step1: modify the environment_lxplus.yml by changing the prefix in environment_lxplus.yml to your local path
  - step2: do 'conda env create -f environment_lxplus.yml' to create the so called "DiHiggs" virtual environment. 
 It would take a while to collect and install all necessary packages
  - step3: check the virtual environment by 'conda env list'. activate the environment by 'conda activate DiHiggs' and disable it by 'conda deactivate'
