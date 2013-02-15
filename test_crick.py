import unittest
import sys
import os
import cPickle
from crick.builders.protein_protein_edges.generic_builder import *
from crick.builders.protein_protein_edges.BioGRID_builder import *
from crick.builders.protein_dna_edges.MotifMap_builder import *
from crick.builders.protein_dna_edges.generic_builder import *
from crick.builders.protein_dna_edges.all_builder import *
from crick.builders.protein_pathway_edges.ProteinPathwayEdge_builder import *
from crick.builders.protein_pathway_edges.generic_builder import *
from crick.const import species
from crick.nodes.pathway import NCIPathway
## object classes
#from crick.objects.DNA import DNA, UCSCGene
#from crick.objects.Protein import UniProtProtein
#from crick.objects.Compound import KEGGMetabolite

## general debugging utils
#from crick.utils.general_utils import wait

class david_list(object):
    def __init__(self):
        self.name='';
        self.genes=[];
class david_collection(object):
    def __init__(self):
        self.samples=[];
        return;
    def load (self):
        for curdirname, subdirnames, curfilenames in os.walk('.'):
            for fname in curfilenames:
                if (fname.split('.')[1])!='txt':
                    continue;
                temp=david_list();
                temp.name=fname.split('.')[0];
                print fname
#get file name
                f=open(fname,'r');
                line=f.readline().lower().split();
                print line
                idx=line.index('to');
#find index at uniport column
                for line in f:
                    temp.genes.append(line.split()[idx]);
                f.close();
                self.samples.append(temp);
        return
    def printinfo(self):
        for sample in self.samples:
            print sample.name, len(sample.genes), sample.genes[0];
        return
class crick_tester(object):
    def __init__(self):
        self.n=network.CrickNetwork(species=species.mm9);
        self.name='';
        self.path='/home/yul13/tmp/'
    def load(self, varin):
        """load from a david list"""
        self.n.add_proteins(varin.genes);
        self.name=varin.name+'_network'
        print self.name, self.n;
        self.pickle(self.name+'loaded')
    def pickle(self,name="default"):
        fout=open(self.path+name+".pkl",'wb');
        print "trying to pickle: ", fout;
        cPickle.dump(self,fout);
        fout.close();
        return
    def unpickle(self,name='/home/yul13/tmp/default.pkl'):
        """unpickles. for some reason manual reasignment is needed, also uses full path"""
        fin=open(name,'rb');
        print "trying to unpickle: ", fin;
        self=cPickle.load(fin);
        print "self now contains", self.n
        fin.close();
        return self
    def exportfig(self,name='testing'):
        """export xgmml file"""
        fig="yu/"+name+".xgmml"
        self.n.export_XGMML(fig,networkout=True)
        return;
    def cleanup(self,deep=False):
        """simple cleanup operations. Do deep at the end"""
        self.n.remove_isoform_proteins();
        self.n.remove_nodes_without_data();
        self.n.remove_similar_trembl_proteins();
        if deep:
            self.n.remove_unconnected_nodes();
        return;
    def closed_ppi(self):
        print "trying to build ppi", self.n
        ppibuilder= BioGRIDProteinProteinEdgeBuilder(self.n);
        ppibuilder.build_network(for_nodes='all', closed_network=True,
                cache=False);
        self.cleanup();
        self.pickle()
        self.exportfig(self.name+'_closed_ppi');
# a simple closed ppi network
        return
    def open_ppi(self):
        print "trying to build ppi", self.n
        ppibuilder= BioGRIDProteinProteinEdgeBuilder(self.n);
        ppibuilder.build_network(for_nodes='all', closed_network=False,
                cache=False);
        self.cleanup();
        self.pickle()
        self.exportfig(self.name+'_open_ppi');
# a simple closed ppi network
        return
    def open_dna(self,domain='source',depth=2):
        # add transcription edges
        print "trying to build motifmap", self.n
        mfbuilder=MotifMapEdgeBuilder(self.n);
        mfbuilder.build_network(upstream = 4000,
                downstream=4000,
                bbls=1.5,
                exon=1,
                fdr=0.5,
                max_depth=depth,
                closed_network=False,
                search_domain=domain,
               #cache=False,
            );
        self.cleanup();
        self.pickle()
        self.exportfig(self.name+'_open_dna_'+domain)
        return
    def closed_dna(self,domain='source', depth=5):
        # add transcription edges
        print "trying to build motifmap", self.n
        mfbuilder=MotifMapEdgeBuilder(self.n);
        mfbuilder.build_network(upstream = 20000,
                downstream=20000,
                bbls=0.5,
                exon=1,
                fdr=1,
                max_depth=depth,
                search_domain=domain,
                closed_network=True,
               #cache=False,
            );
        self.cleanup();
        self.pickle()
        self.exportfig(self.name+'_closed_dna_'+domain)
        return
    def closed_pathway(self):
        """build pathway edges"""
        print "trying to build pathway edges"
        pabuilder=NCIProteinPathwayEdgeBuilder(self.n, target=NCIPathway);
        pabuilder.build_network(closed_network=True,for_nodes='all');
        self.cleanup();
        self.pickle();
        self.exportfig(self.name+'_closed_pathway');
        return
def main():
    c=david_collection();
    c.load();
    c.printinfo();
    d=crick_tester();
    #d=d.unpickle()
    d=d.unpickle('/home/yul13/tmp/DP_1_networkloaded.pkl')
   # d.open_ppi();
    d.closed_dna();
    #d.cleanup(True);
    print d.n
    d.exportfig();
main();

