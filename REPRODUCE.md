# Reproduction instruction for all plots and tables


## Golden set generation and data preparation

1. Run 1.SF1.ipynb;
2. Run 2.DataPreparation.ipynb;
3. Run F1/Figure1.ipynb;


## Proof-of-Concept study (Figure 2)

1. Run ../scripts/F2bc.sh ../data/interim/fig2/;
2. Run ../scripts/F2ade.sh:
    1. ../scripts/F2ade.sh 18 active 10 ../data/interim/fig2/;
    2. ../scripts/F2ade.sh 18 random 10 ../data/interim/fig2/;
3. Run ../scripts/F2a.sh;
4. Run ../scripts/F2f.sh;
    1. ../scripts/F2f.sh 18 2 10 active ../data/interim/fig2/;
    2. ../scripts/F2f.sh 18 3 10 active ../data/interim/fig2/;
    3. ../scripts/F2f.sh 18 4 10 active ../data/interim/fig2/;
    4. ../scripts/F2f.sh 18 5 10 active ../data/interim/fig2/;
    5. ../scripts/F2f.sh 18 6 10 active ../data/interim/fig2/;
    6. ../scripts/F2f.sh 18 7 10 active ../data/interim/fig2/;
    7. ../scripts/F2f.sh 18 8 10 active ../data/interim/fig2/;
    8. ../scripts/F2f.sh 18 9 10 active ../data/interim/fig2/;
5. Run ../scripts/F2_gather_results.sh;
6. Run 3.ResultPreparation_F2bc.ipynb;
7. Run 4.ResultPreparation_F2adef.ipynb;
8. Run 5.F2_T1_SF2_SF3_SF4.ipynb;
9. Run 6.SN1;


## Experimental study (Figure 3 and 4, Supplementary Figures 2, 3, 4, 5, Table 1)

For this part, we've run swiper_cli.py by hand for all experimental xlsx files. The following list is the list of runs, in chronological order: 

1. python ../src/swiper_cli.py -d ../data/external/2024_08_27_SpCas9_set_1_results.xlsx -p ../data/processed/GS_7nt.fasta -s 2024 -q 18 -m RF -o "SpCas9_first.zip"
2. python ../src/swiper_cli.py -d ../data/external/2024_08_29_PpCas9_set_1_results.xlsx -p ../data/processed/GS_7nt.fasta -s 2024 -q 18 -m RF -o "PpCas9_first.zip"
3. python ../src/swiper_cli.py -d ../data/external/2024_10_22_PpCas9_set_2_results.xlsx -p ../data/processed/GS_7nt.fasta -s 2024 -q 18 -m RF -o "PpCas9_second.zip"
4. python ../src/swiper_cli.py -d ../data/external/2024_10_21_SpCas9_set_2_results.xlsx -p ../data/processed/GS_7nt.fasta -s 2024 -q 18 -m RF -o "SpCas9_second.zip"
5. python ../src/swiper_cli.py -d ../data/external/2024_11_11_PpCas9_set_3_results.xlsx -p ../data/processed/GS_7nt.fasta -s 2024 -q 18 -m RF -o "PpCas9_third.zip"
6. python ../src/swiper_cli.py -d ../data/external/2024_11_07_SpCas9_set_3_results.xlsx -p ../data/processed/GS_7nt.fasta -s 2024 -q 18 -m RF -o "SpCas9_third.zip"
7. python ../src/swiper_cli.py -d ../data/external/2024_11_29_PpCas9_set_4_results.xlsx -p ../data/processed/GS_7nt.fasta -s 2024 -q 18 -m RF -o "PpCas9_fourth.zip"
8. python ../src/swiper_cli.py -d ../data/external/2024_11_21_SpCas9_set_4_results.xlsx -p ../data/processed/GS_7nt.fasta -s 2024 -q 18 -m RF -o "SpCas9_fourth.zip"
9. python ../src/swiper_cli.py -d ../data/external/2024_12_06_SpCas9_set_5_results.xlsx -s 2024 -q 18 -m RF -o "SpCas9_fifth.zip" -p None
10. python ../src/swiper_cli.py -d ../data/external/2024_12_09_PpCas9_set_5_results.xlsx -s 2024 -q 18 -m RF -o "PpCas9_fifth.zip" -p None
11. python ../src/swiper_cli.py -d ../data/external/2025_02_19_SpCas9_set_6_results.xlsx -s 2024 -q 18 -m RF -o "SpCas9_sixth.zip" -p None
12. python ../src/swiper_cli.py -d ../data/external/2025_02_21_PpCas9_set_6_results.xlsx -s 2024 -q 18 -m RF -o "PpCas9_sixth.zip" -p None
13. python ../src/swiper_cli.py -d ../data/external/2025_03_18_AlkaCas9_set_1_results.xlsx -s 2024 -q 18 -m RF -o "AlkaCas9_first.zip" -p ../data/processed/GS_7nt.fasta
14. python ../src/swiper_cli.py -d ../data/external/2025_03_24_AlkaCas9_set_2_results.xlsx -s 2024 -q 18 -m RF -o "AlkaCas9_second.zip" -p ../data/processed/GS_7nt.fasta
15. python ../src/swiper_cli.py -d ../data/external/2025_03_25_AlkaCas9_set_3_results.xlsx -s 2024 -q 18 -m RF -o "AlkaCas9_third.zip" -p ../data/processed/GS_7nt.fasta
16. python ../src/swiper_cli.py -d ../data/external/2025_03_27_AlkaCas9_set_4_results.xlsx -s 2024 -q 18 -m RF -o "AlkaCas9_fourth.zip" -p ../data/processed/GS_7nt.fasta
17. python ../src/swiper_cli.py -d ../data/external/2025_04_03_AlkaCas9_set_5_results.xlsx -s 2024 -q 18 -m RF -o "AlkaCas9_fifth.zip" -p None
18. python ../src/swiper_cli.py -d ../data/external/2025_06_05_AlkaCas9_set_6_results.xlsx -s 2024 -q 18 -m RF -o "AlkaCas9_sixth.zip" -p None
19. python ../src/swiper_cli.py -d ../data/external/2025_05_30_SpCas9_set_7_result.xlsx -s 2024 -q 18 -m RF -o "SpCas9_seventh.zip" -p None
20. python ../src/swiper_cli.py -d ../data/external/2025_06_04_PpCas9_set_7_result.xlsx -s 2024 -q 18 -m RF -o "PpCas9_seventh.zip" -p None
21. python ../src/swiper_cli.py -d ../data/external/2025_07_29_AlkaCas9_set_7_results.xlsx -s 2024 -q 18 -m RF -o "AlkaCas9_seventh.zip" -p None
22. python ../src/swiper_cli.py -d ../data/external/2025_07_31_SpCas9_set_8_result.xlsx -s 2024 -q 18 -m RF -o "SpCas9_eigth.zip" -p None
23. python ../src/swiper_cli.py -d ../data/external/2025_08_01_PpCas9_set_8_result.xlsx -s 2024 -q 18 -m RF -o "PpCas9_eigth.zip" -p None
24. python ../src/swiper_cli.py -d ../data/external/2025_09_03_AlkaCas9_set_8_results.xlsx -s 2024 -q 18 -m RF -o "AlkaCas9_eigth.zip" -p None
25. python ../src/swiper_cli.py -d ../data/external/2025_09_08_SpCas9_set_9_results.xlsx -s 2024 -q 18 -m RF -o "SpCas9_ninth.zip" -p None
26. python ../src/swiper_cli.py -d ../data/external/2025_09_10_PpCas9_set_9_results.xlsx -s 2024 -q 18 -m RF -o "PpCas9_ninth.zip" -p None
27. python ../src/swiper_cli.py -d ../data/external/2025_10_27_AlkaCas9_set_9_results.xlsx -s 2024 -q 18 -m RF -o "AlkaCas9_ninth.zip" -p None
28. python ../src/swiper_cli.py -d ../data/external/2025_11_10_SpCas9_set_10_results.xlsx -s 2024 -q 18 -m RF -o "SpCas9_final.zip" -p None
29. python ../src/swiper_cli.py -d ../data/external/2025_11_12_PpCas9_set_10_results.xlsx -s 2024 -q 18 -m RF -o "PpCas9_final.zip" -p None
30. python ../src/swiper_cli.py -d ../data/external/2025_11_25_AlkaCas9_set_10_results.xlsx -s 2024 -q 18 -m RF -o "AlkaCas9_final.zip" -p None

Then unpack each archive in the corresponding directory PROTEIN_NAME/RUN_NUMBER (e.g. AlkaCas9/3), but 10-th run should be renamed to "final". After that, run 7.F3_F4_SF5.ipynb.