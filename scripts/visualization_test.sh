#!/bin/bash


python ../src/swiper_cli.py -d ../data/external/2024_08_27_SpCas9_set_1_results.xlsx -s 23 -q 18 -a -o "test_1" -p ../data/processed/GS_7nt.fasta -t /tmp/swiper_workdir

python ../src/swiper_cli.py -d ../data/external/2024_10_21_SpCas9_set_2_results.xlsx -s 42 -q 18 -a -o "test_2" -p ../data/processed/GS_7nt.fasta -t /tmp/swiper_workdir

python ../src/swiper_cli.py -d ../data/external/2024_11_07_SpCas9_set_3_results.xlsx -s 35 -q 18 -a -o "test_3" -p ../data/processed/GS_7nt.fasta -t /tmp/swiper_workdir

python ../src/swiper_cli.py -d ../data/external/2024_11_21_SpCas9_set_4_results.xlsx -s 55 -q 18 -a -o "test_4" -p ../data/processed/GS_7nt.fasta -t /tmp/swiper_workdir

python ../src/swiper_cli.py -d ../data/external/2024_12_06_SpCas9_set_5_results.xlsx -s 55 -q 18 -a -o "test_5" -p None -t /tmp/swiper_workdir

python ../src/swiper_vis.py --name VisTest --output VisTest --experiments test_1.zip,test_2.zip,test_3.zip,test_4.zip,test_5.zip