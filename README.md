# CompRanking

CompRanking is a pipeline for comprehensively assessing the AMR key message of metagenomic samples and rank their risk of resistome. CompRanking calculate each AMR key message, antibiotic resistance genes abundance, mobility and potentials to acquired by pathogens at the contigs level. Generally, CompRanking can give results derived from three features as the co-occurrance of ARGs, MGEs on one contigs and their potentials in pathogens.

# Getting Started
## Installing
**Step 1:** Change the current working directory to the location where you want the cloned CompRanking directory to be made. **Step 2:** Clone the repository using git command
```sh
git clone https://github.com/GaoyangLuo/CompRanking
```

## Environment settings
### Create environment
Please firstly set up all the environment by the following commands. These commands will help to config all the environment needed.
```sh
$ cd CompRanking
$ conda env create -f CompRanking.yaml
$ bash setup.sh
```
### Setting conda path
CompRanking relies on multi conda environments. Before run the demo test, conda bin path should be pre-requisit. Please set your **absolute bin path** of miniconda. For example, your absolute bin path is `/home/username/miniconda3/bin`.

How to set your absolute conda bin path:

**First step**, vim your `test_yaml.yaml` file
```sh
$ vi test_yaml.yaml
```
**Second step**, Re-write the real path of `miniconda/bin`
```yaml
CompRanking:
  abs_path_to_conda_bin: /your_real_path/miniconda/bin #don't use "~" or "./", please use absolute path
```
**Pleast note that don't use relative path, do not use "~" or "./"**

## Databse download
You can download the databases from the location: `https://doi.org/10.5281/zenodo.8073486`. Or run the command lines below.
```sh
$ wget https://zenodo.org/record/8073486/files/CompRanking_database_v1.tar.gz?download=1
$ wget https://zenodo.org/record/8073486/files/localDB.zip?download=1
$ tar -zxvf CompRanking_database_v1.tar.gz && mv CompRanking_database_v1.tar.gz databases
$ unzip localDB.zip
```

## Demo test
We provided a set of data for test.
```sh
$ python cpr_multiprocess.py -i test_data -t 12 -r 1 -p test_demo
``` 

## Run gene prediction
**Step 1:** Gene prediction can generate contextural information of AMR and pathogen information of the whole metagenome. Run the command line below:
```sh
$ python cpr_multiprocess.py -i <input_dir> -t <threads> -r <if_restart> -p <project_name_prefix>
``` 

Parameters: 
- -i, <input_dir> contains all the fastq files and fasta files. Files of the the sample should be named using identical `<prefix>`. For example, `FileNameOne_1.fq`, `FileNameOne_2.fq` and `FileNameOne.fa` represents the pair-end reads fastq files (after quality control) and the assembly file (containing contigs and pleast do not cut into your customed length, default_min_length=500, which cannot be altered).
- -t, <threads>, the threads you want to use to run the process (Default=16).
- -r, <if_restart> 0 or 1. 0 means continue to run after the last break up point. 1 means re-start from the begeining.
- -p, <project_name_prefix> You should

---

**Step 2:** After finishing all the prediction steps, we should calculate the relative abundance of functional genes, run the command line below:
```sh
$ python ./compranking/GeneCal.py -i <input_dir> -p <project_name_prefix>
```

---

**Step 3:** Generate a risk score and corresponding valuse of each sample. In this step, you can acquire various parameters such as how many ARGs-carried contigs or phage- or plasmids-related contigs in your samples. Please run the command line below:
```sh
$ python ./compranking/baseInfoExtra.py -i <input_dir> -p <project_name_prefix>
```
