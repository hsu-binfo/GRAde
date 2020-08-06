# GRAde
A GRA detection pipeline

## Clone GRAde via Git Large File Storage (LFS) 
Install Git Large File Storage: https://git-lfs.github.com/
```
sudo apt-get install git-lfs 
git clone https://github.com/yihsuanhung/GRAde.git
```

## Install GRAde via miniconda
If you have miniconda or anaconda installed, you can skip this step.
<br>
Install miniconda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html 
<br>
<br>
Install GRAde dependencies with the following command
```
cd GRAde
conda env create -n GRAde --file environment.yml
conda activate GRAde
```
Done!
<br>
## Test your installation
```
python GRAdetection.py 
python GRAdetection.py --fastq example/example1.fastq.gz --outdir /path/to/output/data/ --samplename GRAtesting
```
<br><br>
July 2020
