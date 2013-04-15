#Manages data from cyber-t result (as in typical tables) and out to other managers
#do simple set analyses on top of cyber-t.
import copy
import os
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import sqlite3 as sql
import subprocess
import crick.utils.parsers.data_table as dtable
DBPATH="/home/yul13/data/databases/cybert.db"
LABELS="(Sample TEXT, ProbeID INT, GeneSym TEXT, GeneDescription TEXT,GeneID TEXT, ReferenceID TEXT,\
        Mean REAL, SD REAL, BayesSD REAL, pValDifferential REAL, Bonferonni REAL, BH REAL, DF INT,\
        BayesDF INT, UID INT, RawINFO BLOB)";
BASEDIR="/home/yul13/data/template/"
CURDIR=os.getcwd();
from cybert_slaves import *
VERBOSE=0
ARGS={"FNAME":"","CONTROL":"-1","EXPERIMENT":"-1", "BACKGROUND":"8", "WINDOWSIZE":"101",\
        "INTC":"T","INTE":"T","DOMULT":"T","PPDE":"F","DOQC":"F"}

class cybert_manager(object):
    """cybert t run and result manager,
    manages multiple pairwise comparison with ease,
    currently not for a single multi conditional comparison"""
    def __init__ (self, basedir=CURDIR):
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
            if "order.list" in i:
                bisect=False;
                self.load_ord(i);
                continue;
        if run:
            print "Found raw data, attempting to get sample structures"
            self.loadf();
            if bisect:
                print "No specific order found, assuming Control1-N, Exp1-N"
                self.gen_ord();
        else:
            print "Raw data not found, abort!"
            raise IOError;
        return
    def pickle(self):
        f=open(self.name+".cpkl","wb")
        pickle.dump(self, f);
        return;

    def gen_ord(self):
        self.RUNNAMES=[self.rawname.strip(".xls")]
        self.RUNS=[1];
        self.ORDER=[["C0" for i in xrange(self.numsam/2)]+["E0" for i in xrange(self.numsam/2)]]
        print self.RUNNAMES, self.ORDER, self.RUNS
        return

    def load_ord(self,name):
        f=open("order.list","r")
        i=0;
        self.RUNNAMES=[];self.ORDER=[];self.RUNS=[];
        for line in f:
            if i%2==0:
                self.RUNNAMES.append(line.strip());
            else:
                self.ORDER.append(line.strip().split())
                self.RUNS.append(int(line.strip().split()[-1].strip("E").strip("C")))
            i+=1;

        #now read this order for pairwise comparison
        print self.RUNNAMES, self.ORDER, self.RUNS
        return;

    def loadf(self):
        f=open(self.rawname,"r")
        self.RAWDATA=[];
        for line in f:
            header= [i.strip().strip(".CEL") for i in line.split("\t") if i!=""];
            break;
        print header
        for line in f:
            if line.strip()=="":
                continue;
            temp=line.strip("\n").split("\t");
            assert(len(temp)==len(header)+1)
            self.RAWDATA.append(temp)
        self.numsam=len(header)
        assert(self.numsam%2==0);
        # our preprocessing script does not have column name for probeid
        self.sams=header;
        return;

    def subfile(self,run_num):
        """Make individual CyberT files."""
        controls=[i for i in xrange(len(self.order)) if self.order[i]=="C"+str(run_num)]
        exps=[i for i in xrange(len(self.order)) if self.order[i]=="E"+str(run_num)]
        g=open(str(run_num)+".csv","w")
        g.write("\t".join(["ProbeID"]+[
            self.sams[i] for i in controls+exps
            ]))
        g.write("\n")
        for temp in self.RAWDATA:
            c=[temp[i+1] for i in controls]
            e=[temp[i+1] for i in exps]
            g.write("\t".join([temp[0]]+c+e))
            g.write("\n")
        g.close()
        return;

    def write_runs(self):
        for i in xrange(len(self.RUNNAMES)):
            RUNNAME=self.RUNNAMES[i]
            self.order=self.ORDER[i]
            self.runs=self.RUNS[i]
            try:
                os.mkdir(RUNNAME)
            except OSError:
                pass;
            os.chdir(RUNNAME)
            #now assemble args to insert into template
            template=open(BASEDIR+"cybert_template.R","r").read();
            jobtemplate=open(BASEDIR+"submit.sge","r").read();
            bashtemp="#!/bin/bash \n"
            for i in xrange(self.runs+1):
                myargs=ARGS;
                myargs["FNAME"]=str(i)+".csv"
                self.subfile(i)
                a=self.order.count("C"+str(i))
                b=self.order.count("E"+str(i))
                myargs["CONTROL"]=str(a)
                myargs["EXPERIMENT"]=str(b)
                ARGS["BACKGROUND"]=8-(min([a,b]))
                f=open(str(i)+".R","w"); f.write(template.format(**myargs))
                f.close();
                f=open(str(i)+".sge","w");
                f.write(jobtemplate.format(str(i)+".R",self.run_name))
                f.close();
                bashtemp+="qsub ./{0}\n".format(str(i)+".sge");
                if self.runs==0:
                    break;
            f=open("jobs.sh","w")
            f.write(bashtemp);
            f.close();
            os.chdir(CURDIR)
        return;
    def write_output(self):
        """given output and idconversion file, rejoin the table into more readable form"""
        try:
            os.mkdir("./CyberT_Output")
        except:
            pass;
        f=open("./idConversions.csv","r")
        lookup={};
        for line in f:
            break;
        for line in f:
            ids=line.strip().split("\t")
            try:
                lookup[ids[0]]
                continue;
            except KeyError:
                lookup[ids[0]]=[ids[1],ids[2],ids[3]]
            finally:
                pass;
        for i in xrange(len(self.RUNNAMES)):
            RUNNAME=self.RUNNAMES[i]
            self.order=self.ORDER[i]
            self.runs=self.RUNS[i]
            os.chdir(RUNNAME)
            #these are assumed valid
            for i in xrange(self.runs+1):
                #load output file
                myprefix=str(i)+".csv"
                f=open(myprefix+"cybert.out","r").readlines()
                g=open(CURDIR+"/CyberT_Output/"+self.run_name+"_"+RUNNAME+str(i)+".txt","w")
                g.write("\t".join(["probe_id","gene_sym","gene_name","uniprot"]+f[0].split("\t")))
                controls=[q for q in xrange(len(self.order)) if self.order[q]=="C"+str(i)]
                cnt=0;
                for temp,line  in zip(self.RAWDATA,f[1:]):
                    #do a simple QC
                    temp2=line.split("\t")
                    try:
                        a= "Found row with probe {0}, gene sym: {1} ".format(temp[0],lookup[temp[0]])
                        if VERBOSE: print a;
                        cnt+=1;
                    except KeyError:
                        if VERBOSE:
                            print "Rows with no gene sym ignored {0}".format(temp[0])
                        continue;
                    assert (temp[controls[0]+1]==temp2[0])
                    g.write("\t".join([temp[0]]+lookup[temp[0]]+temp2))
                print "Done merging {0} with {1:d} rows".format(myprefix, cnt)
                g.close()
                print "purging large files with prefix",myprefix
                os.remove(myprefix)
                os.remove(myprefix+"cybert.out")

            os.chdir(CURDIR)

    def make_collections(self,pfx="CD271"):
        """load the output file and make them into cybert_collection objects for later use"""
        pass;
    def serialize(self):
        """serialize all cyber-t collections into databases"""

def main():
    c=cybert_manager()
    c.write_output()


if __name__=="__main__":
    main();


