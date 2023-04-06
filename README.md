# CompRanking: Comprehensively ranking the environmental metagenome samples with resistome assessment

version 1.0

Metagenomic samples resistome assessment and risk ranking pipeline, also known as m3Rs. It's a pipeline for comprehensively ranking the environmental metagenome samples by virtue of assessing their resistome profiles. Generally giving results derived from three features as the co-occurrance of ARGs and MGEs on one contigs and the globle statement of the enrichment of virulence factors.

# workflow
![图片](/image/CompRanking.png)

# How to use
## Environment settings
Please firstly set up all the environment by the following commands. These commands will help to config all the environment needed and a modified seeker command_line.py file.
```sh
$ bash setup.sh
```

## Databse download
```py

```

## Run gene prediction
Step 1:Gene prediction can generate contextural information of AMR and pathogen information of the whole metagenome. Run the command line below:
```sh
$ python cpr_multiprocess.py -i <input_dir> -t <threads> -r <if_restart> -p <project_name_prefix>
``` 

Step 2:After finishing all the prediction steps, we should calculate the relative abundance of functional genes, run the command line below:
```sh
$ python ./compranking/GeneCal.py -i <input_dir> -p <project_name_prefix>
```
Step 3: Generate a risk score and corresponding valuse of each sample. In this step, you can acquire various parameters such as how many ARGs-carried contigs or phage- or plasmids-related contigs in your samples. Please run the command line below:
```sh
$ python ./compranking/baseInfoExtra.py -i <input_dir> -p <project_name_prefix>
```
