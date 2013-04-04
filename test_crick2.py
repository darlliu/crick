import unittest
import sys
import os
import cPickle as pickle
import gzip
import crick
import matplotlib as mpl
mpl.use('Agg')
from crick.builders.protein_protein_edges.generic_builder import *
from crick.builders.protein_protein_edges.BioGRID_builder import *
from crick.builders.protein_dna_edges.MotifMap_builder import *
from crick.builders.protein_dna_edges.generic_builder import *
from crick.builders.protein_dna_edges.all_builder import *
from crick.builders.protein_pathway_edges.ProteinPathwayEdge_builder import *
from crick.builders.protein_pathway_edges.generic_builder import *
from crick.builders.enzyme_metabolite_edges.KEGG_builder import *
from crick.const import species
from read_cybert_table import *
from builder_expander import *
import numpy as np
import matplotlib.pyplot as plt
#from crick.nodes.pathway import NCIPathway
#from crick.viz.google_chart import *
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
        for fname in os.listdir('./'):
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
            print sample.name, len(sample.probes), sample.probes[0],len(sample.genes), sample.genes[0];
        return

class probeset_lookup(object):
    def __init__(self):
        self.lookup={};
        self.db=[];
    def load(self,name="probeset_data.csv"):
        f=open(name,'r');
        #raise exception if file not found
        for line in f:
            if line[0]=='#':continue;
            else:
                temp=line.split('\t');
                for item in temp:
                    self.lookup[item]=temp.index(item);
                break;
        for line in f:
            self.db.append(line.split('\t'));
        f.close();
        print self.lookup;
        return
    def pickle(self):
        f=open("probeset.db",'wb');
        cPickle.dump(self,f);
        f.close();
    def unpickle(self):
        f=open("probeset.db","rb");
        self=cPickle.load(f);
        f.close();
        return self;
    def generate_probe_to_uniprot(self):
        """generate two way lookup table"""
    def generate_swissprot_to_probe(self):
        out={};
        for temp in self.db:
            idfrom=self.lookup["SwissProt"];
            idto=self.lookup["Probe Set ID"];
            keys=temp[idfrom];
            for key in keys.split('///'):
                out[key]=temp[idto];
        return out;

class crick_tester(object):
    def __init__(self):
        self.n=network.CrickNetwork(species=species.mm9);
        self.name='';
        self.path='/home/yul13/tmp/'
        self.descriptions={};
        self.pathwayinfo={};
    def update(self,another):
        assert (type(self)==type(another));
        print "attempting combining", another.name 
        self.n.copy_network(another.n);
        self.name="combined";
        self.path=another.path;
        self.descriptions.update(another.descriptions);
        self.pathwayinfo.update(another.pathwayinfo);
    def combine(self,testers):
        for tester in testers:
            self.update(tester);
        return;
    def load(self, varin):
        """load from a david list"""
        self.name=varin.name+'_network'
        for idx in xrange(len(varin.genes)):
            gene=varin.genes[idx];
            self.n.add_proteins([gene]);
            for item, data in self.n.nodes(data=True):
                try:
                    if data.features["mylookup"]: continue;
                except KeyError:
                    data.features["t_source"]=varin.name.split('_')[0].lower();
                    data.features["mylookup"]=gene;
                    self.descriptions[gene]=varin.descriptions[idx];
                    data.features["t_david_description"]=' '.join(varin.descriptions[idx]);
        self.pickle(self.name+'loaded');
        return;
    def load_refid2(self,varin):
        """load in refids for nodes"""
        for item, data in self.n.nodes(data=True):
            try:
                lookup=data.features["mylookup"];
                try:
                    data.features["probe_refid"]=varin[lookup];
                    print "added ref: ", varin[lookup]
                except KeyError:
                    continue;
            except KeyError:
                lookup=data.features["accessions"];
                try:
                    for key in lookup:
                        data.features["probe_refid"]=varin[key];
                except KeyError:
                    continue;
        retur5;
    def load_refid(self,varin):
        """load in refids for nodes"""
        for item, data in self.n.nodes(data=True):
            try:
                lookup=data.features["mylookup"];
            except KeyError:
                lookup=item;
            try:
                data.features["probe_refid"]=varin.probes[varin.genes.index(lookup)];
            except ValueError:
                continue;
    def pickle(self,name="default"):
        fout=open(self.path+name+".pkl",'wb');
        gout=open("loaded.order",'wa');
        gout.write(name);
        gout.write("\n");
        print "trying to pickle: ", fout;
        cPickle.dump(self,fout);
        fout.close();
        return
    def fix_my_lookup(self):
        for node, data in self.n.nodes(data=True):
            try:
                data.features["mylookup"];
                data.features["mylookup"]=data.features["id"];
                #data.features["t_david_description"]=\
                        #' '.join(self.descriptions[data.features['id']]);
            except KeyError:
                pass;
        return;
    def fix_david_description(self):
        for node, data in self.n.nodes(data=True):
            try:
                data.features["t_david_description"];
            except KeyError:
                data.features["t_david_description"]="";

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
        m=builder_expander();
        ppibuilder= BioGRIDProteinProteinEdgeBuilder(self.n);
        m.modify_builder(ppibuilder);
        ppibuilder.build_network(for_nodes='all', closed_network=True,
                cache=False);
        self.cleanup();
        self.pickle(self.name+"_closedppi")
        self.exportfig(self.name+'_closedppi');
# a simple closed ppi network
        return
    def metabolite(self):
        print "trying to build metabolites"
        kegbuilder=KEGGEnzymeMetaboliteEdgeBuilder(self.n);
        kegbuilder.build_network();
        return;
    def open_ppi(self):
        print "trying to build ppi", self.n
        ppibuilder= BioGRIDProteinProteinEdgeBuilder(self.n);
        ppibuilder.build_network(for_nodes='all', closed_network=False,
                cache=False);
        self.cleanup();
        self.annotate_other_nodes();
        self.pickle(self.name+"_openppi")
        self.exportfig(self.name+'_openppi');
# a simple closed ppi network
        return
    def open_ppi_ones(self):
        print "trying to build ppi", self.n
        ppibuilder= BioGRIDProteinProteinEdgeBuilder(self.n);
        for key,item in self.n.nodes(data=True):
            try:
                if item.features["t_source"]!="other":
                    print "handling core node", key
                    ppibuilder.build_network(for_nodes='all', closed_network=False,
                            cache=False);
                else:
                    pass;
            except KeyError:
                continue;
        self.cleanup();
        self.annotate_other_nodes();
        self.pickle(self.name+"_openppi2")
        self.exportfig(self.name+'_openppi2');
# a simple closed ppi network
        return
    def open_dna_ones(self,domain='target'):
        """for our own nodes, connect with open network"""
        print "trying to build motifmap", self.n
        mfbuilder=MotifMapEdgeBuilder(self.n);
        for key,item in self.n.nodes(data=True):
            try:
                if item.features["t_source"]!="other":
                    print "handling core node", key
                    mfbuilder.build_network(upstream = 3000,
                            downstream=3000,
                            bbls=1,
                            exon=0,
                            fdr=0.5,
                            closed_network=False,
                            search_domain=domain,
                            for_nodes=[key]
                        );
            except KeyError:
                continue;
        self.cleanup();
        self.annotate_other_nodes();
        self.pickle(self.name+"_opendna2")
        self.exportfig(self.name+'_open_dna2_'+domain)
        return
    def open_dna(self,domain='source'):
        # add transcription edges
        print "trying to build motifmap", self.n
        mfbuilder=MotifMapEdgeBuilder(self.n);
        mfbuilder.build_network(upstream = 3000,
                downstream=3000,
                bbls=2,
                exon=0,
                fdr=0.4,
                closed_network=False,
                search_domain=domain,
            );
        self.cleanup();
        self.annotate_other_nodes();
        self.pickle(self.name+"_opendna")
        self.exportfig(self.name+'_open_dna_'+domain)
        return
    def clear_connections(self):
        """method for clearing previous connections annotation"""
        for key , data in self.n.nodes(data=True):
            data.features['_connections']=0
            data.features["_connected_nodes"]=[];
        return;

    def annotate_connections(self):
        """go through each node and annotate its connections, inflating its size at the same time,
        also, for connected nodes count the number of nodes that connects to them,
        as of now we assume directed edges such that no redundant edges are produced"""
        print "annotating connections"
        self.clear_connections();
        for keys in self.n.edge.keys():
            if type(keys)!=type(['']):
                keys=[keys];
            for key in keys:
                print key
                if len(self.n.edge[key])==0: continue;
                else:
                    pair_keys=self.n.edge[key].keys();
                    if type(pair_keys)!=type(['']):
                        pair_keys=[pair_keys];
                try:
                    self.n.node[key].features['_connections']+=len(pair_keys);
                except KeyError:
                    self.n.node[key].features['_connections']=len(pair_keys);
                    print "Key error found at _connections",key
                #annotate self
                for pair_key in pair_keys:
                    try:
                        self.n.node[pair_key].features['_connections']+=1;
                    except KeyError:
                        self.n.node[pair_key].features['_connections']=1;
                    try:
                        self.n.node[key].features['_connected_nodes'].append(pair_key);

                    except KeyError:
                        print "Key error found at _connected_nodes", key
                        self.n.node[key].features['_connected_nodes']=[pair_key];
                    try:
                        self.n.node[pair_key].features['_connected_nodes'].append(key);
                    except KeyError:
                        self.n.node[pair_key].features['_connected_nodes']=[key];
                    if self.n.node[key].features['t_source']!="other":
                        self.n.node[pair_key].features['_connected_to_core']=True;
                    if self.n.node[pair_key].features['t_source']!="other":
                        self.n.node[key].features['_connected_to_core']=True;
        for key, data in self.n.nodes(data=True):
            try:
                data.features['_connections'];
            except KeyError:
                data.features['_connections']=0;
            data.features['_t_inflate']=50+1*data.features['_connections'];
            if data.features['_t_inflate']>150:
                data.features['_t_inflate']=150;
        return;
    def propagate_pathway_info(self):
        """propagate the annotation of pathway information"""
        def do_propagate(ID):
            print "propagating",ID
            try:
                children=self.n.node[ID].features['connected_nodes'];
                self.n.node[ID].features['_noreturn']=True;
            except KeyError:
                return ID;
            for child in children:
                print "now at child", child
                try:
                    self.n.node[child].features['tPathwayInfo']|=self.n.node[ID].features['tPathwayInfo'];
                except KeyError:
                    return ID;
                if self.n.node[child].features['_noreturn']==True:
                    return ID;
                else:
                    return do_propagate(child);
        for key,data in self.n.nodes(data=True):
            data.features['_noreturn']=False;
        for key,data in self.n.nodes(data=True):
            if len(data.features['tPathwayInfo'])==1 and data.features['tPathwayInfo']!=['none']:
                do_propagate(key);
        return

    def annotate_to_keep(self,strict=True,depth_limit=1):
        """annotate whether or not a node is to be kept, the rules are :
            1. if marked t_source not other, is_ff_ or has core pathway info then keep.
            2. if connected to two or more nodes then keep
            3. otherwised marked for delete"""
        for key, data in self.n.nodes(data=True):
            flag=0;
            try:
                if data.features['t_is_ff_from_list']=='false':
                    if data.features['_connections']<=depth_limit:
                        if strict:
                            if data.features['t_source']=='other' \
                                    and data.features["_connected_to_core"]==False:
                                flag=1;
                            else:
                                continue;
                        else:
                            flag=1;
            except KeyError:
                print "There is an error at pruning node ",key,"deleting it anyway"
                flag=1;
            if flag:
                data.features['_tobedeleted']=True;
        return;
    def annotate_prune_added(self,key):
        """Prune nodes that are added except one"""
        for i, j in self.n.nodes(data=True):
            if j.features["t_source"]=="other" and i != key:
                j.features["_tobedeleted"]=True;
            else:
                j.features["_tobedeleted"]=False;
        self.prune();
        return;

    def prune(self):
        """go through nodes delete those marked to be deleted while keeping their keys"""
        try:
            self.deleted;
        except:
            self.deleted=[];
        for key, data in self.n.nodes(data=True):
            try:
                print "attempting delete routine", key
                if data.features["_tobedeleted"]==True:
                    print "deleted unneeded node", key
                    self.n.remove_node(key);
                    self.deleted.append(key);
            except KeyError:
                continue;
        return;
    def closed_dna(self,domain='source'):
        # add transcription edges
        print "trying to build motifmap", self.n
        m=builder_expander();
        mfbuilder=MotifMapEdgeBuilder(self.n);
        m.modify_builder(mfbuilder);
        mfbuilder.build_network(upstream = 4000,
                downstream=4000,
                bbls=1.0,
                exon=0,
                fdr=0.5,
                search_domain=domain,
                closed_network=True,
            );
        self.cleanup();
        self.pickle(self.name+"_closeddna")
        self.exportfig(self.name+'_closed_dna_'+domain)
        return

    def closed_pathway(self):
        """build pathway edges"""
        print "trying to build pathway edges"
        pabuilder=NCIProteinPathwayEdgeBuilder(self.n, target=NCIPathway);
        pabuilder.build_network(closed_network=True,for_nodes='all');
        self.cleanup();
        self.pickle(self.name+"closed_pathway");
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
    

    def annotate_pathway_from_source (self):
        for node, data in self.n.nodes(data=True):
            try:
                annote=[data.features["t_source"]];
                data.features["tPathwayInfo"]=set(annote);
            except KeyError:
                data.features["tPathwayInfo"]=set(['Added']);
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
                        print "found a feature ", lookup, data.features["name"].lower();
                        annote.append(key);
                        break;
            if annote:
                data.features["tPathwayInfo"]=set(annote);
            else:
                data.features["tPathwayInfo"]=set(['none']);

            #print "going through", data.features["name"]
            #if data.features["tPathwayInfo"]:
            #    print "annotation:", data.features["tPathwayInfo"];
        return;
    
    def draw_pathway_chart(self):
        """Simple chart drawing method to annoate pathway source for each node"""
        #template="http://chart.googleapis.com/chart?cht=p&chs=400x200&chd=t:{0}&chl={1}&chtt=Pathway%20Involved&chco={2}"
        template="http://chart.googleapis.com/chart?cht=pc&chtt=Source%20info&chs=200x200&chd=s:C,{0}&chdl={1}&chco={2}|FFFFFF,{3}|FFFFFF"
        template2="http://chart.googleapis.com/chart?cht=pc&chs=200x200&chd=s:C,{0}&chco={1}|FFFFFF,{2}|FFFFFF"
        labels=self.pathwayinfo.keys();
        #labels2=['dp','df','ors','mx','ml','other']
        #labels2=['oe-ko+uniprot','oe+ko-uniprot','phosphatases','phosphatase TF','immune','immune TF']
        #labels2=['phosphatases','phosphatase TF','immune','immune TF']
        labels2=['OE+','OE-','KO+','KO-','OE+KO-','OE-KO+'];
        colors=['004444','FF0000','00FF00','0000FF','FF00FF','FFFF00','00FFFF','880088','008888'];
        colort=['882222','884422','882244','882266','880022','005050'];

        if len(colors)<len(labels):
            raise IndexError;
        mycolordict={};
        mycolordict2={};
        for keyword in labels:
            mycolordict[keyword]=colors[labels.index(keyword)];
        for keyword in labels2:
            mycolordict2[keyword]=colors[labels2.index(keyword)];
        for node, data in self.n.nodes(data=True):
            chartlabels=[];
            chartcolors=[];
            #first we get the inner circle which is the tissue type
            try:
                temp=data.features['t_source'];
                chartlabels.append(temp);
                chartcolors.append(mycolordict2[temp]);
            except KeyError:
                chartlabels.append('unknown');
                chartcolors.append('EEEEEE');
            for pathway in data.features["tPathwayInfo"]:
                try:
                    chartcolors.append(mycolordict[pathway]);
                    chartlabels.append('%20'.join(pathway.split()));
                except KeyError:
                    continue;
            numbers=['C' for i in xrange(len(chartlabels)-1)]; 
            if len(chartlabels)==1:
                data.features["tChartUrl"]=self.wrap(template.format("2",chartlabels[0]+"|N/A",chartcolors[0],"FFFFFF"));
                data.features["_CustomGraphics"]=template2.format("2",chartcolors[0],"FFFFFF");
                continue;
            else:
                data.features["tChartUrl"]=self.wrap(template.format(\
                        "".join(numbers), '|'.join(chartlabels),chartcolors[0],'|'.join(chartcolors[1:])\
                        ));
                data.features["_CustomGraphics"]=template2.format(\
                        "".join(numbers), chartcolors[0],'|'.join(chartcolors[1:])\
                        );
                print data.features["tChartUrl"];
        return
    
    def annotate_other_nodes(self):
        for nodes,data in self.n.nodes(data=True):
            try:
                data.features["t_source"]
                continue;
            except KeyError:
                data.features["t_source"]="other";
        return;

    
    def annotate_tf_from_list(self):
        """check if the gene is in a given gene list"""
        
        return;
    #def annotate_cybert_results_hair(self):
        #f=gzip.open("means.gz","rb");
        #means=pickle.load(f);
        #g=gzip.open("bsds.gz","rb");
        #bsds=pickle.load(g);
        #h=gzip.open("sds.gz","rb");
        #sds=pickle.load(h);
        ##h=open("labels.lb","rb");
        ##labels=pickle.load(h);
        #labels=['MX','ORS','DF','DP','ML']
        #for node, data in self.n.nodes(data=True):
            #try:
                #mean=means[data.features["probe_refid"]];
                #bsd=bsds[data.features["probe_refid"]];
                #sd=sds[data.features["probe_refid"]];
##                print "got results from cybert", len(mean), len(bsd)
                #data.features["cybert_means"]=mean;
                #am=self.draw_hist(mean,bsd,sd,labels);
                #data.features["cybert_sd_bayes"]=bsd;
                #data.features["cybert_plot"]=am;
            #except:
                #data.features["cybert_means"]=[];
                #data.features["cybert_sd_bayes"]=[];
        #return
    #def merge_tfbs_network(self,fname="./pipeline-TF.pkl"):
        #f=open(fname,"rb");
        #a=pickle.load(f);
        #lookup={};
        #for key, data in self.n.nodes(data=True):
            #try:
                #lookup[key]=data.features["t_source"]
            #except KeyError:
                #data.features["t_source"]="other";
        #self.n.copy_network(a);
        #for key, source in lookup.items():
            #self.n.node[key].features["t_source"]=source;
        #return;

    def annotate_cybert_results(self):
        C1=cybert_collection();
        C2=cybert_collection();
        C1.loaddb("./cybert_OEvsWT1cybert.db");
        C2.loaddb("./cybert_KOvsWT2cybert.db");
        #labels=["WT1","OE", "WT2","KO"]
        labels=["WT1","OE", "WT2","KO"]
        #labels=self.labels;
        for node, data in self.n.nodes(data=True):
            try:
                c1wt=C1.samples[0].entry_by_genesym(data.features["name"]);
                c1oe=C1.samples[1].entry_by_genesym(data.features["name"]);
                c2wt=C2.samples[0].entry_by_genesym(data.features["name"]);
                c2ko=C2.samples[1].entry_by_genesym(data.features["name"]);
                sname="";
                if c1wt.mean>c1oe.mean and c1wt.pval<0.05:
                    data.features["tPathwayInfo"]|=set(["OE_DOWN"]);
                    sname+="OE-"
                elif c1wt.mean<c1oe.mean and c1wt.pval<0.05:
                    data.features["tPathwayInfo"]|=set(["OE_UP"]);
                    sname+="OE+"
                if c2wt.mean>c2ko.mean and c2wt.pval<0.05:
                    data.features["tPathwayInfo"]|=set(["KO_DOWN"]);
                    sname+="KO-"
                elif c2wt.mean<c2ko.mean and c2wt.pval<0.05:
                    data.features["tPathwayInfo"]|=set(["KO_UP"]);
                    sname+="KO+"
                if sname:
                    data.features["t_source"]=sname;
                if type(c1wt)==type(None):
                    raise KeyError;
                print c1wt,data.features["name"]
                mean=[c1wt.mean, c1oe.mean,c2wt.mean,c2ko.mean];
                bsd=[c1wt.bsd+mean[0], c1oe.bsd+mean[1],c2wt.bsd+mean[2],c2ko.bsd+mean[3]];
                sd=[c1wt.sd+bsd[0], c1oe.sd+bsd[1],c2wt.sd+bsd[2],c2ko.sd+bsd[3]];
#                print "got results from cybert", len(mean), len(bsd)
                data.features["cybert_means"]=mean;
                data.features["cybert_pvals"]=[c1wt.pval,c2wt.pval];
                data.features["IlluminaProbe"]=c1wt.refid;
                am=self.draw_hist(mean,bsd,sd,labels);
                data.features["cybert_sd_bayes"]=bsd;
                data.features["cybert_sd"]=sd;
                data.features["cybert_plot"]=am;
                print "annotated one cybert result", node
            except:
                print "error getting cybert info at key", node
                data.features["cybert_means"]=[];
                data.features["cybert_sd_bayes"]=[];
                data.features["cybert_pvals"]=[];
                data.features["cybert_sd"]=[];
        return
    def draw_hist(self,means,bsds,sds,labels,name="CyberT%20Means%20and%20Bayesian%20SD"):
        """draw a double layer histogram given data"""
        template="https://chart.googleapis.com/chart?cht=bvo&chs=300x200&chtt={4}&chd=t:{0}|{1}|{2}&chdl=bsd|sd|mean&chco=4D89F9,C6D9FD,4DD622&chxt=x,y&chxl=0:|{3}&chbh=a&&chds=a&chxr=1,{5:.3f}{6:.3f}"
        start=min(means)*0.98;
        end=max(sds)*1.02;

        r1=[str(entry-start) for entry in bsds];
        r11=[str(entry-start) for entry in sds];
        r2=[str(entry-start) for entry in means];
        return self.wrap(template.format(",".join(r1),",".join(r11),",".join(r2),'|'.join(labels),name,start,end));
    
    def wrap(self,url):
        """wraps url in img link"""
        return """<img src="{0}"></img>""".format(url);

    def annotate_tf_from_descriptions(self,keywords=\
            ["transcription","core","factor","binding","regulat","polymerase",\
             "DNA"]):
        """attemp to annotate tf identity from david information"""
        for key,item in self.n.nodes(data=True):
            for lookup in keywords:
                try:
                    if lookup in self.descriptions[key] \
                            or lookup in " ".join(item.features["comment_similarity"]):
                        item.features["t_is_ff_from_list"]="true";
                        item.features["border_color"]="#FF6600"
                        break;
                except KeyError:
                    item.features["t_is_ff_from_list"]="false";
        for key,item in self.n.nodes(data=True):
            try:
                item.features["t_is_ff_from_list"];
            except KeyError:
                item.features["t_is_ff_from_list"]="false";
            try:
                if item.features["border_color"]=="#ffff00":
                    item.features["t_is_ff_from_list"]="true";
            except KeyError:
                continue;
        return;
    def do_annotation_routine(self):
        self.get_pathway_source();
        self.annotate_other_nodes();
        self.annotate_pathway();
        self.annotate_tf_from_descriptions();
        self.annotate_connections();
        self.annotate_to_keep();
        self.prune();
        self.propagate_pathway_info();
        #self.annotate_cybert_results();
        self.draw_pathway_chart();
    def load_david_clusters(self,fname):
        """generate david annotation clusters from output file"""
        try:
            self.clusters;
        except:
            self.clusters={};
        clusters=[];
        temp={};
        f=open(fname);
        flag=0;
        enrich=0.0;
        for line in f:
            if line.strip()=="":
                clusters.append(temp);
                print "added cluster of length {0:d}".format(len(temp));
                flag=0;
                temp={};
                enrich=0.0;
                continue;
            elif "Annotation Cluster" in line:
                lookup={};
                enrich=float(line.split()[-1])
                print "found cluster enrichment of {0:.3f}".format(enrich);
                flag=1;
                continue;
            elif flag==1:
                header=[i.lower().strip() for i in line.split("\t") if i!=""];
                for i in header:
                    lookup[i]=header.index(i);
                flag=2;
                continue;
            elif flag==2:
                items=line.split("\t");
                goterms=items[lookup["term"]];
                go=goterms.split('~')[0];
                meaning=goterms.split('~')[-1];
                temp[go]=[enrich, meaning,items[lookup["genes"]],items];
                continue;
            else:
                raise IOError;

        self.clusters[fname]=clusters;
        print "added {0:d} clusters".format(len(clusters));
        return
    def generate_david_pie_charts(self):
        """generate a very simple pie chart from david annotation clustering"""
        noadd=["positive","process","and","or","the","negative","pathway","regulates","function"]
        def wordpair(wordlists):
            freq={};
            for i in wordlists:
                try:
                    freq[i]+=1;
                except KeyError:
                    freq[i]=0;
            keys=freq.keys();
            vals=freq.values();
            first=keys[vals.index(max(vals))]
            vals.pop(vals.index(max(vals)))
            keys.pop(keys.index(first))
            second=keys[vals.index(max(vals))]
            return [first,second]
        assert len(self.clusters)>0;
        for clustername,clusters in self.clusters.items():
            plt.figure(figsize=(10,10));
            plt.axes([0.1,0.1,0.8,0.8])
            enrich=[];
            annotes=[];
            for cluster in clusters:
                meanings=[];
                for i, j in cluster.items():
                    meanings+=[i for i in j[1].split() if len(i)>2 and i not in noadd];
                en=j[0];
                keywords="\n".join(wordpair(meanings))
                print "enrichment {0:.3f} of {1}".format(en,keywords);
                enrich.append(en);
                annotes.append(keywords);
            print enrich
            if len(enrich)>10:
                enrich=enrich[:10]
                annotes=annotes[:10]
            plt.pie(enrich, labels=annotes);
            plt.savefig(clustername+".png")
        return;
    def export_uniprots(self):
        f=open(self.name+"_uniprots.txt","w")
        for key in self.n.nodes():
            f.write(key);
            f.write("\n");
        return;
    def export_probe(self):
        f=open(self.name+"_probes.txt","w")
        for key,data  in self.n.nodes(data=True):
            try:
                f.write(data.features["IlluminaProbe"].upper());
            except KeyError:
                f.write(key);
            f.write("\n");
        return;
    def lookup_enrichment_from_tfbs_pipeline(self):
        """look up and add tfbs enrichment info from the pipeline"""
        try:
            fname=self.FNAME;
        except:
            fname="default.txt"
        f=open("./tfbs/"+fname,"r");
        lookup={};
        #for line in f:
            #header=[i.strip("\"") for i in line.split() if i!=""];
            #break;
        #for i in header:
            #lookup[i]=header.index(i);
        #print lookup;
        self.tfbs={};
        for line in f:
            temp=[i.strip("\"").upper() for i in line.split() if i!=""];
            print temp
            #print temp[lookup["motif.uniprots"]],\
                    #temp[lookup["motif.accessions"]],temp[lookup["report.odds"]]
            self.tfbs["".join(temp[0].split("-")).upper()]=[temp[5],\
                    temp[6],temp[7],temp[9],\
                    temp];
            self.tfbs[temp[5]]=[temp[5],\
                    temp[6],temp[7],temp[9],\
                    temp];
        for key,item in self.n.nodes(data=True):
            try:
                temp=self.tfbs[key];
                print "found a tfbs enrichment at ", key
                item.features["TFBS_EnrichScore"]=str(temp[2]);
                item.features["TFBS_EnrichPval"]=str(temp[3]);
                item.features["TFBS_Acce"]=temp[1];
            except KeyError:
                try:
                    temp=self.tfbs[item.features["name"]];
                    print "found a tfbs enrichment at ", key
                    item.features["TFBS_EnrichScore"]=str(temp[2]);
                    item.features["TFBS_EnrichPval"]=str(temp[3]);
                    item.features["TFBS_Acce"]=str(temp[1]);
                except KeyError:
                    item.features["TFBS_EnrichScore"]="";
                    item.features["TFBS_EnrichPval"]="";
                    item.features["TFBS_Acce"]="N/A"
        return;
    def lookup_tissue_source_from_david(self):
        f=open("./txts/david_tissue.txt")
        for line in f:
            header=[i.lower().strip() for i in line.split("\t") if i.strip()!=""];
            break;
        self.tissue={};
        for line in f:
            temp=[i.upper().strip() for i in line.split("\t") if i.strip()!=""];
            prots=[i.strip() for i in temp[header.index("id")].split(",")];
            for i in prots:
                self.tissue[i]=temp[3];
        for key,item in self.n.nodes(data=True):
            try:
                ii=item.features["IlluminaProbe"];
            except:
                item.features["tissue"]="not significant"
                continue;
            if ii in self.tissue.keys():
                item.features["tissue"]=self.tissue[ii];
            else:
                item.features["tissue"]="not significant"
    def annotate_prune_except(self, key):
        exempt=self.n.node[key].features["_connected_nodes"];
        for i, j in self.n.nodes(data=True):
            if i ==key:
                continue;
            if i not in exempt:
                j.features["_tobedeleted"]=True;
            else:
                j.features["_tobedeleted"]=False;
        self.prune();
        return;

    def write_pretty_pipeline_report(self):
        return
    def export_summary(self):
        """for each node in the network, export a table of its genesym,
        source, whether or not a TF, enrichment of motifs, david description/
        comment when available and then its connection in the network."""
        template="\t".join(["GeneSym","Uniprot","Source of addition",\
                "tissue","comment1","comment2","is TF","TFBS Enrichment",\
                "Binding Motif accession","network connections",\
                "WT1 mean","OE mean","WT2 mean","KO mean",
                "WT1 BayesSD","OE BayesSD","WT2 BayesSD","KO BayesSD",
                "WT1 SD","OE SD","WT2 SD","KO SD",
                "OE pval", "KO pval","Pathway Info"
                ]);
        f=open("summary_"+self.name+".xls","w");
        f.write(template)
        f.write("\n");
        for key, item in self.n.nodes(data=True):
            means=[str(i) for i in item.features["cybert_means"]];
            bsds=[str(i) for i in item.features["cybert_sd_bayes"]];
            sds=[str(i) for i in item.features["cybert_sd"]];
            pvals=[str(i) for i in item.features["cybert_pvals"]];
            print means, bsds,sds,pvals
            try:
                line="\t".join([\
                item.features["name"],key,item.features["t_source"],\
                item.features["tissue"],item.features["t_david_description"],\
                ",".join(item.features["comment_similarity"]),\
                item.features["t_is_ff_from_list"],item.features["TFBS_EnrichScore"],\
                item.features["TFBS_Acce"],str(item.features["_connections"])]+\
                means+bsds+sds+pvals+[",".join(list(item.features["tPathwayInfo"]))]\
                )
                f.write(line);
                f.write("\n");
            except KeyError:
                print "an error occured at ", key
                continue;
        return;
    def do_export_routine(self):
        self.annotate_cybert_results();
        self.fix_david_description();
        self.annotate_tf_from_descriptions();
        self.lookup_tissue_source_from_david();
        self.lookup_enrichment_from_tfbs_pipeline();
        self.export_summary();
        return;

def load_all_crick_objects(basedir='/home/yul13/tmp/load/'):
    cricks=[];
    for name in os.listdir(basedir):
        temp=crick_tester();
        temp=temp.unpickle(basedir+name);
        cricks.append(temp);
    return cricks;
def initialize_all_crick_objects():
    c=david_collection();
    c.load();
    c.printinfo();
    cricks=[];
    for sample in c.samples:
        d=crick_tester();
        d.load(sample);
        cricks.append(d);
    return cricks;
def load_refids(cricks):
    c=david_collection();
    c.load();
    c.printinfo();
    e=probeset_lookup();
    e.load();
    f=e.generate_swissprot_to_probe();
    for d in cricks:
        d.load_refid2(f);
        for sample in c.samples:
            d.load_refid(sample);
    return cricks;

def main():
    networks=load_all_crick_objects();
    #networks=initialize_all_crick_objects();
    for i in networks:
        i.closed_dna();
        i.closed_ppi();
        i.exportfig(i.name+'initial');
    return
    #networks=load_refids(networks);
    #combined=crick_tester();
    #combined.combine(networks);
    #combined.do_annotation_routine();
    #combined.exportfig(combined.name+"combined_temp");
    #combined.pickle(combined.name+"combined_temp")
    #d=d.unpickle('/home/yul13/tmp/ORS_2_networkloaded.pkl')
    #d=d.unpickle()
    #d.annotate_pathway();
    #d.annotate_tf_from_descriptions();
    #d.draw_pathway_chart();
#    print d.pathwayinfo;
    #d.open_ppi();
    #d.closed_dna();
    #d.cleanup(True);
    #print d.n
#main()
