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
from crick.viz.google_chart import *
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
        self.pathwayinfo={};
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
                    data.features["t_source"]=varin.name.split('_')[0].lower();
                    data.features["mylookup"]=gene;
                    self.descriptions[gene]=varin.descriptions[idx];
                    data.features["t_david_description"]=' '.join(varin.descriptions[idx]);
        self.pickle(self.name+'loaded')
    def pickle(self,name="default"):
        fout=open(self.path+name+".pkl",'wb');
        gout=open("loaded.order",'wa');
        gout.write(name);
        gout.write("\n");
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
        self.annotate_other_nodes();
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
        self.annotate_other_nodes();
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
        f=open("coregenes.list",'r');
        self.pathwayinfo={};
        temp=[];
        for line in f:
            line=line.strip();
            if line=='': continue;
            if '#' in line : continue;
            if line[0]=='%':
                self.pathwayinfo[line[1:].strip()]=temp;
                temp=[];
                continue;
            temp.append(line.strip().lower());
        return
    def annotate_pathway (self):
        """a very simple method for annotating\
                which pathway the genes come from
                """
        for node, data in self.n.nodes(data=True):
            annote=[];
            for key,item in self.pathwayinfo.items():
                for lookup in item:
                    if lookup in data.features["name"].lower():
                        annote.append(key);
                        break;
            data.features["tPathwayInfo"]=annote;
            #print "going through", data.features["name"]
            #if data.features["tPathwayInfo"]:
            #    print "annotation:", data.features["tPathwayInfo"];
        return;
    def draw_pathway_chart(self):
        """Simple chart drawing method to annoate pathway source for each node"""
        #template="http://chart.googleapis.com/chart?cht=p&chs=400x200&chd=t:{0}&chl={1}&chtt=Pathway%20Involved&chco={2}"
        template="http://chart.googleapis.com/chart?cht=p&chs=400x200&chd=t:{0}&chl={1}&chco={2}"
        labels=self.pathwayinfo.keys();
        colors=['000000','FF0000','00FF00','0000FF','FF00FF','FFFF00','00FFFF','888888','222222'];
        if len(colors)<len(labels):
            raise IndexError;
        mycolordict={};
        for keyword in labels:
            mycolordict[keyword]=colors[labels.index(keyword)];
        for node, data in self.n.nodes(data=True):
            chartlabels=[];
            chartcolors=[];
            for pathway in data.features["tPathwayInfo"]:
                try:
                    chartcolors.append(mycolordict[pathway]);
                    chartlabels.append('%20'.join(pathway.split()));
                except KeyError:
                    continue;
            numbers=['2' for i in xrange(len(chartlabels))]; 
            if len(chartlabels)==0:
                data.features["tChartUrl"]=template.format("2","None","000000");
                continue;
            else:
                data.features["tChartUrl"]=template.format(\
                        ",".join(numbers), '|'.join(chartlabels),'|'.join(chartcolors)\
                        );
                data.features["_CustomGraphics"]=data.features["tChartUrl"];
                print data.features["tChartUrl"];

        return
    def annotate_other_nodes(self):
        for nodes,data in self.n.nodes(data=True):
            try data.features["t_source"]:
                continue;
            except KeyError:
                data.features["t_source"]="other";
        return;

    def annotate_tf_from_list(self):
        """check if the gene is in a given gene list"""
        
        return;
    def annotate_tf_from_descriptions(self,keywords=\
            ["transcription","core","factor","binding","regulat","polymerase"\
                    ]):
        """attemp to annotate tf identity from david information"""
        for key,item in self.n.nodes(data=True):
            for lookup in keywords:
                if lookup in self.descriptions[key]:
                    item.features["t_is_ff_from_list"]="true";
                    break;
                else:
                    item.features["t_is_ff_from_list"]="false";
def main():
    c=david_collection();
    c.load();
#    c.printinfo();
    d=crick_tester();
    for sample in c.samples:
        d.load(sample);
        d.annotate_pathway();
        d.annotate_tf_from_descriptions();
        d.draw_pathway_chart();
    #d=d.unpickle()
    d=d.unpickle('/home/yul13/tmp/ORS_2_networkloaded.pkl')
    d.get_pathway_source();
#    print d.pathwayinfo;
    #d.open_ppi();
    #d.closed_dna();
    #d.cleanup(True);
    #print d.n
    d.exportfig("mapping");
main();

