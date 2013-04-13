#Manages data from cyber-t result (as in typical tables) and out to other managers
#do simple set analyses on top of cyber-t.
import copy
import os
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import sqlite3 as sql
import crick.utils.parsers.data_table as dtable
DBPATH="/home/yul13/data/databases/cybert.db"
LABELS="(Sample TEXT, ProbeID INT, GeneSym TEXT, GeneDescription TEXT,GeneID TEXT, ReferenceID TEXT,\
        Mean REAL, SD REAL, BayesSD REAL, pValDifferential REAL, Bonferonni REAL, BH REAL, DF INT,\
        BayesDF INT, UID INT, RawINFO BLOB)";
BASEDIR=os.getcwd();
from cybert_slaves import *

class cybert_manager(object):
    """cybert t run and result manager"""
    def __init__ (self, basedir=BASEDIR):
        self.col=cybert_collection();
        self.run_name=basedir.split('/')[-1];
        tempfiles=os.listdir("./")
        run=False
        bisect=True;
        for i in tempfiles:
            if "processed_data" in i:
                run=True;
                self.rawname=i;
                continue
            if "conditions" in i:
                bisect=False;
                self.load_ord(i);
                continue;
        if run:
            print "Found raw data, attempting to get sample structures"
            if bisect:
                print "No specific order found, assuming Control1-N, Exp1-N"
                self.gen_ord();
        else:
            print "Raw data not found, abort!"
            raise IOError;





