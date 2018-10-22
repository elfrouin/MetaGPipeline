#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 10:51:49 2016

@author: frouin.e
"""
"""
Calculates expression per ko per sample. Outputs a tsv with ko expression in
each sample. There are 2 columns. The header looks as follows

ko_hit\thist_mean_sum_namesample

Parses output from
    - GHOSTKOALA output file for ko functional annotation on proteins.aa.fa
    - histogram files for each sample generated with
        prodigal -i assembly.fa -f gff -p meta > prodigal.gff
        bedtools coverage -hist samplen.bam prodigal.gff > samplen.gff.coverage.hist
"""
import argparse
import csv
from itertools import groupby


class Histo():
    def __init__(self,line):
        description = line.split("\t")
        self.ID = description[8].split(";")[0][3:]
        self.hist_cov = int(description[9])
        self.hist_ratio = float(description[12])


class KoAnnot():
    def __init__(self,line):
        annotation = line.strip("\n").split("\t")
        id_name = annotation[0].split("contig-100_")[1].split("_")
        self.ID = str(int(id_name[0])+1)+"_"+id_name[1]
        self.K = annotation[1]


def main(file_k, file_hist, sample_name, file_out):
    # create a list of Histo objects
    hist = []
    with open(file_hist) as file:
        for l in file :
            if not l.startswith('all'):
                hist.append(Histo(l))

    # create a dict with key= ID of the contig and value= coverage of this contig
    hist_dict = {}
    for key, group in groupby(hist, key=lambda x: x.ID):
        hist_dict[key] = sum([obj.hist_cov*obj.hist_ratio for obj in (list(group))])

    # create a list of K_annot objects
    ko = []
    with open(file_k) as file:
        for l in file :
                ko.append(KoAnnot(l))

    # create a dict with key= Ko and value= contigs with this Ko annotation
    ko_dict = {}
    ko.sort(key=lambda x: x.K, reverse=True)
    for key, group in groupby(ko, key=lambda x: x.K):
        ko_dict[key] = [obj.ID for obj in (list(group))]

    # create a dict with key= Ko and value= sum of coverage of contigs with this Ko annotation
    cov_dict = {}
    for key, IDs in ko_dict.items():
        list_cov=[]
        for ID in IDs:
            if ID in hist_dict:
                list_cov.append(hist_dict[ID])
        cov_dict[key] = sum(list_cov)

    # write the final dict in a csv file
    with open(file_out, 'w') as csv_file:
        writer = csv.writer(csv_file, delimiter='\t')
        writer.writerow(["ko_hit", "hist_mean_sum_"+sample_name])
        for key in sorted(cov_dict):
            writer.writerow([key, round(cov_dict[key], 7)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=
                                     argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sampleK", help="GHOSTKOALA output file for ko "
                        "annotation of proteinsfa")
    parser.add_argument("--sampleHist", help="bedtools "
                        "coverage generated file one for each sample")
    parser.add_argument("--out", help="Output file (.tsv)")
    parser.add_argument("--sampleName", default='sample', help="File with sample "
                        "names, one line each. Should be same nr as"
                        "sample.coverage.hist.")
    args = parser.parse_args()

    main(args.sampleK, args.sampleHist, args.sampleName,args.out)
