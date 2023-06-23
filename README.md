# CompRanking
## Comprehensively ranking the environmental metagenome samples with resistome assessment

CompRanking is a pipeline for comprehensively assessing the AMR key message of metagenomic samples and rank their risk of resistome. CompRanking calculate each AMR key message, antibiotic resistance genes abundance, mobility and potentials to acquired by pathogens at the contigs level. Generally, CompRanking can give results derived from three features as the co-occurrance of ARGs, MGEs on one contigs and their potentials in pathogens.

# Getting Started
## Installing
**Step 1:** Change the current working directory to the location where you want the cloned CompRanking directory to be made. **Step 2:** Clone the repository using git command
```sh
git clone https://github.com/GaoyangLuo/CompRanking
```

## Environment settings
Please firstly set up all the environment by the following commands. These commands will help to config all the environment needed.
```sh
$ cd CompRanking
$ conda env create -f CompRanking.yaml
$ bash setup.sh
```

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

**Step 2:** After finishing all the prediction steps, we should calculate the relative abundance of functional genes, run the command line below:
```sh
$ python ./compranking/GeneCal.py -i <input_dir> -p <project_name_prefix>
```
**Step 3:** Generate a risk score and corresponding valuse of each sample. In this step, you can acquire various parameters such as how many ARGs-carried contigs or phage- or plasmids-related contigs in your samples. Please run the command line below:
```sh
$ python ./compranking/baseInfoExtra.py -i <input_dir> -p <project_name_prefix>
```
