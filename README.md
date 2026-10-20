# Swiper

<img src="swiper.png" alt="logo" width="150">


## Sequencing-free PAM discovery for CRISPR-Cas9 orthologs via active learning.

---
Swiper identifies the protospacer-adjacent motif (PAM) of an uncharacterized Cas9 ortholog from a small number of iterative cleavage assays guided by Active Machine Learning.

The Monte-Carlo statistical test shows when the consensus is stable and the search can be completed: Swiper, no swiping!

---
### Quick start

To install, run the following:
```bash

git clone https://github.com/bakirillov/swiper.git
cd swiper
pip install -r requirements.txt
```
After that create an excel file with list named 1 and the following structure as on the example below:
| # | PAM | Sequence | Outcome |
|---|-----|----------|---------|
| 0 | 3   | CTTCCGC  | -       |
| 1 | 4   | ATATGGG  | -       |
| 2 | 5   | AACGGGG  | -       |
| 3 | 6   | TAAGTGT  | -       |

Replace the sequences with your sequences and in ''Outcome'' input the results of your cleavages, ''-'' if no cleavage happened and ''+'' if cleavage is successful.

After that, run:
```bash

python ./src/swiper_cli.py -d YOURFILE.xlsx -p "./supplementary/Supplementary Data 3.fasta" -s YOUR_SEED -q YOUR_ROUND_SIZE -a -o "OUTPUT_MODEL_ARCHIVE_1" -t /tmp/swiper_workdir/

#Use -p None if you want to run against all possible k-mers of given k; 
#Use -p "./supplementary Data 3.fasta" if you are looking for 7nt candidate; 
#Use -p "./supplementary Data 2.fasta" if you are looking for 6nt candidate; 
#If you want your own Golden Set, generate it according to ./notebooks/EDF1.ipynb;
```
This will add list '2' to the original file with the next round sequences and save the model to OUTPUT_MODEL_ARCHIVE_1.zip. 

Perform cleavage assay and fill list '2' with the outcomes.  Repeat until convergence.

To visualize the results and identify convergence, run the following:

```bash

python swiper_vis.py --name YOUR_EXPERIMENT_NAME --output REPORT_NAME --experiments OUTPUT_MODEL_ARCHIVE_1.zip,OUTPUT_MODEL_ARCHIVE_2.zip,OUTPUT_MODEL_ARCHIVE_3.zip...and so on
```

This will create file REPORT_NAME.pdf


---
### Reproduction instructions

Read and follow REPRODUCE.md.
