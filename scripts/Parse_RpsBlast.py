#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 10:33:06 2016

@author: frouin.e
"""

# library
import argparse
import os
from Bio.Blast import NCBIXML
import pandas
import operator

 

# define input file 
parser = argparse.ArgumentParser(description= "Parse blast results (xml file) after prodigal queries in rps.blast against COG database")
parser.add_argument("input_file" ,   help="Input file (xml format)")
parser.add_argument("-COGfunction" , help="Output file : table of ORF counts in each COG number")
parser.add_argument("-COGcategory",  help="Output file : table of ORF counts in each COG category")
parser.add_argument("-partialgene",  help="Apply a weighting for partial genes found by prodigal ", type = bool , default = False)
parser.add_argument("-COG_db",  help="path to a list of COG ID ", default = os.environ['HOME']+'/Documents/01-Databases/rpsblast_db/cognames2003-2014.tab')

args = parser.parse_args()


##############################################################################
##                              CLASS COG
##############################################################################  

class COG :
    def __init__ (self):
        self.ID = ''
        self.cat = []
        # with partial gene weighting        
        self.count = 0
        # without partial gene weighting
        self.count_all = 0
    def new_COG(self,dico,key):
        self.ID= key
        # add in a list one or several categories (for one COG)
        self.cat = list(dico[key][0])
    def uncharacterized (self) :
        self.ID='uncharacterized'
        self.cat = ['uncharacterized']
        
##############################################################################
##                         CLASS COGs_inQuery
##############################################################################  

class COGs_inQuery:
    def __init__(self):
        self.collection = []
        tmp_data =pandas.read_csv(args.COG_db, sep='\t', header =0)
        self.dico_COG = tmp_data.set_index('# COG').T.to_dict('list')
        self.count_queries =0
        list_cat= map(chr, range(65, 91))
        list_cat.append("uncharacterized")
        # dico : COG cat : [count (partial), count_all]
        self.count_cat={ i : [0,0] for i in list_cat }
        
    def create_COGs(self):
        # add each COG (object) to the list self.collection        
        for key in self.dico_COG:
            tmp_COG= COG()
            tmp_COG.new_COG(self.dico_COG,key)
            self.collection.append(tmp_COG)
        # add the uncharaterized object at the en of the list self.collection
        unchar_COG= COG()
        unchar_COG.uncharacterized()
        self.collection.append(unchar_COG)
        self.collection.sort(key=operator.attrgetter('ID'))
        


    def update_db_with_queries(self):
        for record in NCBIXML.parse(open(args.input_file)) :
            # We count the submitted queries             
            self.count_queries +=1
            # We want to ignore any queries with no search results:
            c= 0
            if record.alignments :
#                print "QUERY: %s..." % record.query[:180]
                for align in record.alignments :
                        for hsp in align.hsps :
                            if c == 0 : 
                                tmp_COG= align.hit_def.split(',')[0]
#                               print tmp_COG
                                COG_index =next((i for i, COG in enumerate(self.collection) if COG.ID == tmp_COG), None)
                                if COG_index != None:                            
                                    self.collection[COG_index].count_all +=1   
                                else:
                                    print 'Error: ID %s is not in %s' %(tmp_COG, args.COG_db)
                            else :
                                break
                            c+=1 
            else :
                self.collection[-1].count_all +=1
    
    def COGfunction_out(self) :
        f = open(args.COGfunction ,'w')        
        f.write(("COG\thits_%s\n") %(os.path.basename(args.input_file)))
        for obj in self.collection:
            f.write(('%s\t%i \n') %(obj.ID ,obj.count_all))
        f.close()
            
        
    def COGcategory_out(self) :
        for obj in self.collection:
                for c in obj.cat:
                    # MAJ count_all
                    self.count_cat[c][1]+=(obj.count_all/float(len(obj.cat)))
        f = open(args.COGcategory ,'w')        
        f.write(("COG_cat\thits_%s\n") %(os.path.basename(args.input_file)))
        for key in self.count_cat:
            f.write(('%s\t%f \n') %(key ,self.count_cat[key][1]))
        f.close()
        
        
##############################################################################
##                              MAIN
##############################################################################        


DB_COG= COGs_inQuery()
DB_COG.create_COGs()

DB_COG.update_db_with_queries()
print DB_COG.count_queries

DB_COG.COGfunction_out()
DB_COG.COGcategory_out()
#print DB_COG.dico_COG
