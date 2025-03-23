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
$ conda activate CompRanking
$ pip install MicrobeCensus
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
$ tar -zxvf CompRanking_database_v1.tar.gz && mv CompRanking_database_v1 databases
$ unzip localDB.zip
```

## Demo test
We provided a set of data for test.
```sh
$ python cpr_multiprocess.py -i test_data -t 4 -r 1 -p test_demo
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

**Step 1.5: Calculating scg and AGS** Before Step2 we, we need to generate AGS and scg files.

To generate AGS files, following the command lines below
```sh
cd $PATH_TO_CompRanking
vi $PATH_TO_CompRanking/scripts/AGS.sh #change your input_dir in the workdir
bash $PATH_TO_CompRanking/scripts/AGS.sh
```
---

**Step 2:** After finishing all the prediction steps, we should calculate the relative abundance of functional genes, run the command line below:
```sh
$ python ./compranking/multiGeneCal_metagenome_rpkg_scg_geneName.py 
        -i <input_dir> 
        -p <project_prefix> 
        -n AGS
        -t 16
        -d <pth2KK2db> #this option is for cell copy normalized by sequence abundance, need to run multiGeneCal_16s.py
```
The output demo is like below:

| ARG_name | Class   | Database | MGE_type |Sample_1 |Sample_2 |Sample_3 |
|----------|---------|----------|----------|---------|---------|---------|       
| AAC(2')-I  | aminoglycoside  | DeepARG | Unknown       |0.003 |0.004 |0.006 |
| ERMB       | macrolide       | RGI     | phage/plasmid |0.002 |0.003 |0.004 |
| SUL3       | sulfonamide     | SARG    | plasmid       |0.001 |0.003 |0.005 |

**note**: `phage/plasmid` means ARGs both fpundto be co-located with phage- or plasmid-like contig in one sample (microbial community).
`plasmid` means only found to be co-located with plasmid-like contig. `Unknown` means not to be found co-located with any MGEs, but not representing it is not co-located with any MGEs, probably due to the accuracy and recall of identification method.

---


**Step 3:** Generate a risk score and corresponding valuse of each sample. In this step, you can acquire various parameters such as how many ARGs-carried contigs or phage- or plasmids-related contigs in your samples. Please run the command line below:
```sh
$ python ./compranking/baseInfoExtra_nContigs.py -i <input_dir> -p <project_name_prefix>
```

<div style="overflow-x: auto;">
  <table>
    <tr>
      <th>sample_name/index</th>
      <th>nContigs</th>
      <th>nARGs_contigs</th>
      <th>nMGEs_contig</th>
      <th>nMGEs_plasmid_contig</th>
      <th>nMGEs_phage_contigs</th>
      <th>nPAT_contigs</th>
      <th>nARGs_MGEs_contig</th>
      <th>nARGs_MGEs_plasmid_contigs</th>
      <th>nARGs_MGEs_phage_contigs</th>
      <th>nARGs_MGEs_PAT_contigs </th>
      <th>fARG </th>
      <th>fMGE </th>
      <th>fMGE_plasmid </th>
      <th>fMGE_phage </th>
      <th>fPAT</th>
      <th>fARG_MGE </th>
      <th>fARG_MGE_plasmid </th>
      <th>fARG_MGE_phage</th>
      <th>fARG_MGE_PAT</th>
      <th>score_pathogenic</th>
      <th>score_phage</th>
      <th>score_plasmid</th>
      <!-- Add more columns as needed -->
    </tr>
    <tr>
      <td>Sample2</td>
      <td>433650 </td>
      <td> 482</td>
      <td> 327780</td>
      <td> 270746</td>
      <td> 32238</td>
      <td> 397572</td>
      <td> 357</td>
      <td> 304</td>
      <td> 28</td>
      <td> 61</td>
      <td> 0.0011</td>
      <td> 0.7558</td>
      <td> 0.6247</td>
      <td> 0.0746</td>
      <td> 0.9168</td>
      <td> 0.0008</td>
      <td> 0.0007</td>
      <td> 6.4568</td>
      <td> 0.0001</td>
      <td> 23.1492</td>
      <td> 20.7351</td>
      <td> 22.7372</td>
      <!-- Add more data as needed -->
      <tr>
      <td>Sample1</td>
      <td>433650 </td>
      <td> 482</td>
      <td> 327780</td>
      <td> 270746</td>
      <td> 32238</td>
      <td> 397572</td>
      <td> 357</td>
      <td> 304</td>
      <td> 28</td>
      <td> 61</td>
      <td> 0.0011</td>
      <td> 0.7558</td>
      <td> 0.6247</td>
      <td> 0.0746</td>
      <td> 0.9168</td>
      <td> 0.0008</td>
      <td> 0.0007</td>
      <td> 6.4568</td>
      <td> 0.0001</td>
      <td> 23.1492</td>
      <td> 20.7351</td>
      <td> 22.7372</td>
      <!-- Add more rows as needed -->
    </tr>
    <!-- Add more rows as needed -->
  </table>
</div>

## How to calculate each ARG class and their carriers counts

Use the jupyter notebook `MGE_carried_ARGs_type_count.ipynb` to calculate. The metadata record the number of five types of elements that co-exist with ARGs: plasmid, phage, unclassified (can be any type of sequences, including chromosome or other unknown or unidentified MGEs), IS (Insertion Sequence), IE (Integrated Elements). Table will be generated like this:

|         |Sample1_x |Sample2_y |Sample3_z |Sample3_m |Sample3_n |
|---------|----------|----------|----------|----------|----------|       
| aminoglycoside  |33 |35 |36 |37 |38 |
| macrolide       |44 |45 |46 |37 |38 |
| tetracycline    |22 |23 |24 |37 |38 |

```sh
sampleName_x: #plasmid
sampleName_y: #phage
sampleName_z: #unclassified
sampleName_m: #IS
sampleName_n: #IE
```

## Re-running if pipeline halted
Every process will generte a checkpointing file in the repo checkdone, with file name like `<Your_Project_Name>.index_build.done`. If you want to re-run the pipeline from the last broken step, you can set the parameter `-r` as `0`, which means don't re-run from the beginning. If you set `1`, means you want to re-run from the beginning. You can also delete the `.done` file if you want to re-run the speicific step. We make this pipeline able to identify which step you have run and which one is not completed.
