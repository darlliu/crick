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
        self.probes=[];
        self.descriptions=[];
class david_collection(object):
    def __init__(self):
        self.samples=[];
        return;
    def load (self):
        for curdirname, subdirnames, curfilenames in os.walk('.'):
            for fname in curfilenames:
                if (fname.split('.')[-1])!='txt':
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
                    if line.split()[0] in temp.probes:
                        continue
                    temp.probes.append(line.split()[0]);
                    temp.genes.append(line.split()[idx]);
                    try:
                        temp.descriptions.append\
                                ([item.strip().lower() for item in line.split()[4:]])
                    except IndexError:
                        #in case no description
                        pass;
                f.close();
                self.samples.append(temp);
        return
    def printinfo(self):
        for sample in self.samples:
            print sample.name, len(sample.genes), sample.genes[0],sample.descriptions;
        return
class crick_tester(object):
    def __init__(self):
        self.n=network.CrickNetwork(species=species.mm9);
        self.name='';
        self.path='/home/yul13/tmp/'
        self.descriptions={};
    def load(self, varin):
        """load from a david list"""
        self.name=varin.name+'_network'
        for idx in xrange(len(varin.genes)):
            gene=varin.genes[idx];
            self.n.add_proteins([gene]);
            for item, data in self.n.nodes(data=True):
                try:
                    if data["mylookup"]: continue;
                except KeyError:
                    data.features["tSource"]=varin.name.split('_')[0];
                    data.features["mylookup"]=gene;
                    self.descriptions[gene]=varin.descriptions[idx];
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
    def open_dna(self,domain='target'):
        # add transcription edges
        print "trying to build motifmap", self.n
        mfbuilder=MotifMapEdgeBuilder(self.n);
        mfbuilder.build_network(upstream = 2000,
                downstream=2000,
                bbls=2,
                exon=0,
                fdr=0.4,
                closed_network=False,
                search_domain=domain,
            );
        self.cleanup();
        self.pickle()
        self.exportfig(self.name+'_open_dna_'+domain)
        return
    def closed_dna(self,domain='source'):
        # add transcription edges
        print "trying to build motifmap", self.n
        mfbuilder=MotifMapEdgeBuilder(self.n);
        mfbuilder.build_network(upstream = 2000,
                downstream=2000,
                bbls=1.5,
                exon=0,
                fdr=0.5,
                search_domain=domain,
                closed_network=True,
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
    def get_pathway_source(self):
        f=open("pathway.info",'r');
        self.pathwayinfo={};
        for line in f:
            if line[0]=='#': continue;
            if line[0]=='~':
                self.pathwayinfo[line[1:].strip()]=temp;
                temp=[];
            temp.append(line.strip().lower());
        return
    def annotate_pathway (self):
        """a very simple method for annotating\
                which pathway the genes come from
                """
        return;
    def annotate_tf_from_list(self):
        """check if the gene is in a given gene list"""
        
        return;
    def annotate_tf_from_descriptions(self,keywords=\
            ["transcription","factor","binding","regulate","polymerase"\
                    ]):
        """attemp to annotate tf identity from david information"""
def main():
    c=david_collection();
    c.load();
    c.printinfo();
    d=crick_tester();
    #for sample in c.samples:
    #    d.load(sample);
    #d=d.unpickle()
    d=d.unpickle('/home/yul13/tmp/DP_1_networkloaded.pkl')
    
   # d.open_ppi();
    #d.closed_dna();
    #d.cleanup(True);
    #print d.n
    #d.exportfig();
main();

