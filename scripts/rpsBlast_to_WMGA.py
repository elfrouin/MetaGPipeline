#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 11:51:00 2016

@author: frouin.e
"""

# library
import argparse
import os
from Bio.Blast import NCBIXML
import pandas
 

# define input file 
parser = argparse.ArgumentParser(description= "Parse blast results (xml file) after prodigal queries in rps.blast against COG database")
parser.add_argument("input_file" ,   help="Input file (xml format)")
parser.add_argument("-output" , help="Output file : table WebMGA output like",default= 'rpsBlast_parsed.csv')
parser.add_argument("-COG_db",  help="path to a list of COG ID ", default = '/root/databases/COG/listcogs.txt')
parser.add_argument("-COG_categories",  help="path for a file with COG categories description", default = '/root/databases/COG/COG_categories.txt')

args = parser.parse_args()

        
##############################################################################
##                              MAIN
##############################################################################        
COG_db =pandas.read_csv(args.COG_db,
                        sep='\t',
                        header=None,
                        names =['c1','c2','c3','class', 'short description','COG','long description'],
                        usecols=['class', 'short description','COG','long description'])
COG_cat =pandas.read_csv(args.COG_categories, sep='\t', header =0)
dict_COG_cat = COG_cat.set_index('class').T.to_dict('list')
        
count_queries =0
f = open(args.output ,'w')       
# WARNING score = bits score 
f.write('#Query\tHit\tE-value\tIdentity\tScore\tQuery-start\tQuery-end\tHit-start\tHit-end\tHit-length\tdescription\tclass\tclass description\n')
for record in NCBIXML.parse(open(args.input_file)) :
    # We count the submitted queries             
    count_queries +=1
    # We want to ignore any queries with no search results:
    c= 0
    if record.alignments :
#   print "QUERY: %s..." % record.query[:180]
        for align in record.alignments :
            for hsp in align.hsps :
                if c == 0 : 
                    COG= align.hit_def.split(',')[0]
                    COG_class = COG_db[COG_db.COG==COG]["class"].values[0]
                    identity_percent = str(int(round(hsp.identities/float(hsp.align_length),2)*100))
                    if len(COG_class)<=1:
                        class_description = dict_COG_cat[COG_class][0]
                    else :
                        class_description = 'Multiple classes' 
                    str_cog= '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'  %(record.query.split('#')[0].strip(),COG,hsp.expect,identity_percent,hsp.bits,hsp.query_start,hsp.query_end,hsp.sbjct_start,hsp.sbjct_end,hsp.align_length,COG_db[COG_db.COG==COG]["long description"].values[0],COG_class,class_description)
                    f.write(str_cog)

                else :
                    break
                c+=1 
f.close()
