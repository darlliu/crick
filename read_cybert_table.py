#!/home/baldig/shared_libraries/centos64/pkgs/python/2.6.5/bin/python
import copy
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
import sqlite3 as sql
DBPATH="/home/yul13/data/databases/cybert.db"
LABELS="(Sample TEXT, ProbeID INT, GeneSym TEXT, GeneDescription TEXT,GeneID TEXT, ReferenceID TEXT,\
        Mean REAL, SD REAL, BayesSD REAL, pValDifferential REAL, Bonferonni REAL, BH REAL, DF INT,\
        BayesDF INT, UID INT, RawINFO BLOB)";
class filtered_data(object):
    """simple adaptor for holding filtered out data"""
    def __init__ (self):
        self.name="";
        self.gsm='';
        self.data=set([]);
        return;
class filtered_collection (object):
    """collection of filtered data"""
    def __init__(self):
        self.entries=[];
        self.temp=set([]);
        return;
    def load (self,fname="GSE3142_highly_expressed_by_uniref.txt"):
        f1=open(fname,'r');
        temp2=filtered_data();
        for line in f1:
            if len(line.split('!!'))==2:
                header=line.split('!!');
                temp2.name=header[0];
                temp2.gsm=header[1];
            elif line.strip()=='':
                continue;    
            else:
                temp=line.split();
                for parts in temp:
                    temp2.data.add(parts);
                self.entries+=copy.deepcopy([temp2]);
                header=[];
                temp2.data=set([]);
        f1.close();
        return;
    def find_all_pairwise_union (self):
        for a in self.entries:
            for b in self.entries:
                temp=a.data&b.data;
                print a.name, " and ", b.name, " has common highly expressed genes", len(temp)
        return;
    def find_all_pairwise_unique(self):
        for a in self.entries:
            for b in self.entries:
                temp=a.data-b.data;
                print a.name, " minus ", b.name, " has uniquely highly expressed genes", len(temp)
        return;        
    def find_all_unique_genes(self):
        for a in self.entries:
            temp=a.data;
            for b in self.entries:
                if(a.name!=b.name):
                    temp=temp-b.data;
            print a.name, " has uniquely highly expressed genes", len(temp)
        return;
    def savedb(self):
        f=gzip.open("filtered.db","wb");
        pickle.dump(self.entries,f);
        f.close();
        return;
    def loaddb(self,fname="filtered.db"):
        f=gzip.open(fname,"rb");
        self.entries=pickle.load(f);
        return;

class cybert_entry(object):
    def __init__(self):
        self.refid="";
        self.geneid="";
        self.genedecription="";
        self.genesym="";
        self.refseq="";
        self.mean=0;
        self.sd=0;
        self.bsd=0;
        self.uid=0;
        self.bh=0;
        self.bf=0;
        self.df=0;
        self.dfb=0;
        self.pval=0;
        self.raw={};
        #here we used number as uid
    def serialize(self,cur,table,sample_name):
        """serialize into a database given the cursor"""
        nums=','.join([
            str(self.mean),str(self.sd),str(self.bsd),str(self.pval),str(self.bf),str(self.bh),\
            str(self.df),str(self.dfb),str(self.uid)
                ])
        sentence="INSERT INTO {0} VALUES('{7}','{1}','{2}','{3}','{4}','{5}',{6},'{8}')".format(
                    table,self.refid,self.genesym,self.geneid,self.genedecription,self.refseq\
                            ,nums,sample_name,";".join(self.raw.values())
                    );
        cur.execute(sentence);
        return;
    def give_uid(self):
        return self.uid-1;
    def c(self, another):
        if (self.mean*another.mean*self.bsd*another.bsd==0):
            return 0;
        #print self.bsd,another.bsd
        return (self.mean-another.mean)/((self.bsd**2+another.bsd**2))**(0.5);
    def d(self, another,p=0.05):
        assert self==another;
        assert self.pval==another.pval
        if self.pval  > p:
            return 0 ;
        elif self.mean>another.mean:
            if "ptpn" in self.genesym.lower():
                print self.mean, another.mean,self.pval, self.genesym,self.refid
            return 1;
        else:
            return -1;
    def __eq__(self,another):
        return self.refid==another.refid;
    def __hash__(self,s=None):
        if type(s)==type(None):
            s=self;
        ord3 = lambda x : '%.3d' % ord(x)
        return int(''.join(map(ord3, s.refid)))
class cybert_sample(object):
    def __init__(self):
        self.name="";
        self.entries=[];
        self.cutoff=2.0;
        #herer we assume a cutoff for fairly large df
    def serialize(self,cur,table):
        """serialize all entires"""
        for entry in self.entries:
            entry.serialize(cur,table,self.name);
        return;
    def __eq__(self,another):
        return self.name==another.name;
    def give_entries(self):
        for entry in self.entries:
            yield entry;
    def uids (self):
        for entry in self.entries:
            yield entry.uid;
    def refs (self):
        for entry in self.entries:
            yield entry.refid;
    #def uid_by_cutoff(self, another,cutoff=2.0):
        #"""filter entries by cutoff, output is uid"""
        #for entry in self.entries:
            #if self==another:
                #yield entry.uid;
                ##avoid self filtering
            #elif entry.c(another.entries[entry.give_uid()])>cutoff:
                #yield entry.uid;

    def entry_by_cutoff(self, another,cutoff=5e-2):
        """filter entries by cutoff, output is uid"""
        for entry in self.entries:
            if self==another:
                yield entry;
                #avoid self filtering
            elif entry.d(another.entries[entry.give_uid()],cutoff)==1:
                yield entry;
    def entry_by_cutoff_lo(self, another,cutoff=5e-2):
        """filter entries by cutoff, output is uid"""
        for entry in self.entries:
            if self==another:
                yield entry;
                #avoid self filtering
            elif entry.d(another.entries[entry.give_uid()],cutoff)==-1:
                yield entry;
        return

    def entry_by_genesym(self,key):
        for i in self.entries:
            if key.lower()==i.genesym.lower():
                return i;
            else:
                continue;
        return
    def entry_by_genesym2(self,keys):
        out=[];
        for i in self.entries:
            for key in keys:
                if key.lower() in i.genesym.lower():
                    out.append(i);
                else:
                    continue;
        return out;
    def ref_by_cutoff(self, another,cutoff=5e-2):
        """filter entries by cutoff, output is uid"""
        for entry in self.entries:
            if self==another:
                yield entry.refid;
                #avoid self filtering
            elif entry.d(another.entries[entry.give_uid()],cutoff)==1:
                yield entry.refid;
    def ref_by_cutoff_lo(self, another,cutoff=5e-2):
        """filter entries by cutoff, output is uid"""
        for entry in self.entries:
            if self==another:
                yield entry.refid;
                #avoid self filtering
            elif entry.d(another.entries[entry.give_uid()],cutoff)==-1:
                yield entry.refid;
    def entry_by_uid(self, uid):
        """yield entry by uid (number), used to connect to filtered-data"""
        for entry in self.entries:
            if entry.uid==uid:
                return entry;
    def entry_by_refid(self,refid):
        """yield entry for refid"""
        for entry in self.entries:
            if entry.refid.strip()==refid.strip():
                return entry;

class cybert_collection(object):
    def __init__(self):
        self.samples=[];
        self.name="";
        self._samples=None;
    def load_gene_syms(self):
        g1=open("./probeid/IlluminaProbeLookup.pkl","rb");
        DB=pickle.load(g1);
        for sam in self.samples:
            for en in sam.entries:
                try:
                    en.genesym=DB[en.refid.upper()]["Symbol"];
                    en.refseq=DB[en.refid.upper()]["RefSeq_ID"];
                except KeyError:
                    print "a probe was not found to have gene symbol..."
                    en.genesym="";
                    en.refseq="";
        return;
                    
    def load (self,fname):
        self.name=fname.split('/')[-1].split('.')[0];
        f1=open(fname,'r');
        line=f1.readline();
        labels=[i.strip("\"").lower() for i in line.split() if i];
        print labels;
        lookup={};
        num_sam=len([i for i in labels if "mean" in i]);
        if self._samples!=None:
            assert len(self._samples)==num_sam;
        print "Number of samples (pairwise)", num_sam
        if num_sam >2:
            print "use multisample version please"
            raise IndexError;
        for i in xrange(len(labels)):
            lookup[labels[i]]=i;
            ##    gt_udx=i;
            #if item.lower().find("dfbetbayes")!=-1:
                #df_idx=i;
            #elif item.lower().find("dfwithbayes")!=-1:
                #dfb_idx=i;
            #elif item.lower().find("bonferroni")!=-1:
                #bf_idx=i;
            #elif item.lower().find("bh")!=-1:
                #bh_idx=i;
            #elif item.lower().find("mean")!=-1:
                #num_sam+=1;
                ##get number of samples
                #mean_idx.append(i);
            #elif (item.lower().find("bayessd")!=-1):
                #bsd_idx.append(i);
            #elif(item.lower().find("sd")!=-1):
                #sd_idx.append(i);
        mean_idx=[labels.index(i) for i in labels if "mean" in i];
        sd_idx=[labels.index(i) for i in labels if "rasd" in i];
        bsd_idx=[labels.index(i) for i in labels if "bayessd" in i];
        print mean_idx, sd_idx, bsd_idx;
        at=1;
        for i in xrange(num_sam):
            self.samples+=copy.deepcopy([cybert_sample()]);
            if self._samples==None:
                self._samples=["WT","KO"];
            self.samples[i].name=self._samples[i];
        for line in f1:
            toget=[i.strip().lower() for i in line.split() if i!=""];
            temp=cybert_entry();
            for i in xrange(len(labels)):
                temp.raw[labels[i]]=toget[i];
            temp.uid=at;
            temp.refid=toget[lookup["lab_0"]].strip("\"");
            #temp.genedescription=toget[gt_idx];
            temp.bf=float(toget[lookup["bonferroni"]]);
            temp.bh=float(toget[lookup["bh"]]);
            #temp.df=toget["df"];
            temp.dfb=float(toget[lookup["bayesdf"]]);
            temp.pval=float(toget[lookup["pval"]]);
            means=[];
            sds=[];
            bsds=[];
            for i in xrange(num_sam):
                self.samples[i].entries+=copy.deepcopy([temp]);
                try:
                    means.append((float(toget[mean_idx[i]])));
                except ValueError:
                    print "error parsing",toget[mean_idx[i]];
                    means.append(0);
                try:
                    sds.append((float(toget[sd_idx[i]])));
                except ValueError:
                    print "error parsing",toget[sd_idx[i]];
                    sds.append(0);
                try:
                    bsds.append((float(toget[bsd_idx[i]])));
                except ValueError:
                    print "error parsing",toget[bsd_idx[i]];
                    bsds.append(0);
            #print means
            for i in xrange(num_sam): 
                self.samples[i].entries[-1].mean=means[i];
                self.samples[i].entries[-1].sd=sds[i];
                self.samples[i].entries[-1].bsd=bsds[i];
            at+=1;
        self.printinfo();
    def serialize(self):
        """serialize into database"""
        con=sql.connect(DBPATH);
        
        with con:
            cur=con.cursor();
            cur.execute("DROP TABLE IF EXISTS {0}".format(self.name))
            cur.execute("CREATE TABLE IF NOT EXISTS {0}{1}".format(self.name,LABELS))
#create table if not present and delete all rows
            for sample in self.samples:
                sample.serialize(cur,self.name);
        return;
    def printinfo(self):
        for i in self.samples:
            print i.name
            print i.entries[0].mean,i.entries[0].sd,i.entries[0].pval, i.entries[0].refid;
        return;

    def find_differential_entries(self,sample,lo=False,p=0.05):
        """gives differential entries that are different from the input"""
        out=set(sample.give_entries());
        for pair in self.samples:
            if lo:
                out=out&set(sample.entry_by_cutoff_lo(pair,p));
            else:
                out=out&set(sample.entry_by_cutoff(pair,p));
        return out;
    def find_differential_refids(self,sample,lo=False,p=0.05):
        """gives differential entries that are different from the input"""
        out=set(sample.refs());
        for pair in self.samples:
            if lo:
                out=out&set(sample.ref_by_cutoff_lo(pair,p));
            else:
                out=out&set(sample.ref_by_cutoff(pair,p));
        return out;
    def find_all_unique_uids(self):
        f1=open("unique_uids.txt",'w');
        for sample in self.samples:
            line=sample.name;
            temp=self.find_differential_entries(sample);
            for entry in temp:
                line=line+"\t"+str(entry.uid);
            print "sample:",sample.name," has uniques: ", len(temp);
            f1.write(line);
            f1.write("\n");
        f1.close();
        return    
    def find_all_unique_entries(self,lo=False):
        if lo:
            f1=open("diff_under_refs.txt",'w');
        else:
            f1=open("diff_over_refs.txt",'w');
        for sample in self.samples:
            line=sample.name;
            temp=self.find_differential_entries(sample);
            for entry in temp:
                line=line+"\t"+str(entry.refid);
            if lo:
                print "sample:",sample.name," has differential uniques (underexpressed): ", len(temp);
            else:
                print "sample:",sample.name," has differential uniques (overexpressed): ", len(temp);
            f1.write(line);
            f1.write("\n");
        f1.close();

    def find_all_unique_refs(self,lo=False):
        if lo:
            f1=open("diff_under_refs.txt",'w');
        else:
            f1=open("diff_over_refs.txt",'w');
        for sample in self.samples:
            line=sample.name;
            temp=self.find_differential_entries(sample);
            for entry in temp:
                line=line+"\t"+str(entry.refid);
            if lo:
                print "sample:",sample.name," has differential uniques (underexpressed): ", len(temp);
            else:
                print "sample:",sample.name," has differential uniques (overexpressed): ", len(temp);
            f1.write(line);
            f1.write("\n");
        f1.close();
        return       

    def find_all_unique_differential_refs(self,idss,lo=False):
        """find all differential genes that are also in set ids (refids)"""
        if lo:
            f1=open("unique_low_refs.txt",'w');
        else:
            f1=open("unique_high_refs.txt",'w');
        for idx in xrange(len(self.samples)):
            sample=self.samples[idx];
            #this sample
            ids=idss[idx];
            #set of high refids
            selfids=set(sample.refs());
            #all refs in this sample
            uni_ids=selfids&ids;
            #only the highs in this sample probably redundant
            temp2=set([]);
            for refid in uni_ids:
                temp2=temp2|set(sample.entry_by_refid(refid));
                #get the entries by refid not necessary but just in case
            line=sample.name;
            temp=self.find_differential_entries(sample,lo);
            #all of the unique entries, change this line to make it low
            temp=temp&temp2;
            for entry in temp:
                line=line+"\t"+str(entry.refid);
            if lo:
                print "sample:",sample.name," has highly expressed genes that are uniquely underexpressed: ", len(temp);
            else:
                print "sample:",sample.name," has highly expressed genes that are uniquely overexpressed: ", len(temp);
            f1.write(line);
            f1.write("\n");
        f1.close();
        return       
    def find_all_unique_overexpressed_refs_by_uid(self,idss):
        """find all differential genes that are also in set ids (uids)"""
        f1=open("unique_high_refs.txt",'w');
        for idx in xrange(len(self.samples)):
            sample=self.samples[idx];
            ids=idss[idx];
            selfids=set(sample.uids());
            uni_ids=selfids&ids;
            temp2=set([]);
            for number in uni_ids:
                temp2=temp2|set(sample.entry_by_uid(number));
                #get union of highly expressed genes
            line=sample.name;
            temp=self.find_unique_entries(sample);
            temp=temp&temp2;
            for entry in temp:
                line=line+"\t"+str(entry.refid);
            print "sample:",sample.name," has highly expressed uniques: ", len(temp);
            f1.write(line);
            f1.write("\n");
        f1.close();
        return       
    def savedb(self):
        print self.name
        f=open(self.name.split(".")[0]+"cybert.db","wb");
        pickle.dump(self.samples,f);
        f.close();
        return;
    def loaddb(self,fname="cybert.db"):
        f=open(fname,"rb");
        self.samples=pickle.load(f);
        return self.samples;
    def mergedb(self, fname, samples=None):
        """merge a preloaded db to do more work"""
        f=gzip.open(fname,"rb");
        temp=pickle.load(f);
        if samples!=None:
            assert len(temp)==len(samples);
            for i in xrange(len(temp)):
                temp[i].name=samples[i];
        if self.samples!=[]:
            assert len(self.samples[0].entries)!=len(temp[0].entries)
        self.samples+=temp;
        return temp;
    def export_all_cybert_mean_bsds(self):
        """this is just a method to export a complete lookup table for all
        expression data"""
        means={};
        bsds={};
        sds={};
        for ref in self.samples[0].refs():
            collection=[sample.entry_by_refid(ref) for sample in self.samples];
            if len (collection)!=len(self.samples): raise IndexError;
            temp_mean=[entry.mean for entry in collection];
            temp_bsd=[entry.bsd for entry in collection];
            temp_sd=[entry.sd for entry in collection];
            means[ref]=temp_mean;
            bsds[ref]=temp_bsd;
            sds[ref]=temp_sd;
            #print temp_mean,temp_bsd;
        f=gzip.open("means.gz","wb");
        g=gzip.open("bsds.gz","wb");
        h=gzip.open("sds.gz","wb");
        pickle.dump(means,f);
        pickle.dump(bsds,g);
        pickle.dump(sds,h);
        f.close();
        g.close();
        h.close();
    def plotpvals(self):
        p=[];
        ppde=[];
        cppde=[];
        for i in self.samples[0].entries:
            p.append(i.pval);
            ppde.append(float(i.raw["ppde.p"]))
            cppde.append(float(i.raw["cum.ppde.p"]))
        plt.figure();
        plt.subplot(221);
        plt.hist(np.array(p),bins=30);
        plt.title("p values histogram")
        plt.subplot(222);
        plt.hist(np.array(ppde),bins=30);
        plt.title("ppde.p values histogram")
        plt.subplot(223);
        plt.hist(np.array(cppde),bins=30);
        plt.title("cum.ppde.p values histogram")
        plt.savefig("./hist.png");

def sirt1():
    """special analysis for sirt1"""
    s=cybert_collection();
    p=1e-2;
    s.loaddb("./cybert_WT1vsWT2cybert.db");
#load the WT comparison as background
    assert len(s.samples)==2;
    s.samples[0].name="WT1(OE)"
    s.samples[1].name="WT2(WT)"
    print s.samples[0].entries[0].genesym;
    wt1_up=s.find_differential_entries(s.samples[0],p=p);
    wt1_down=s.find_differential_entries(s.samples[0],True,p=p);
    temp=s.find_differential_entries(s.samples[1],p=p);
    print "number of entries that are up in WT1 is", len (wt1_up)
    print "number of entries that are up in WT2 is", len (wt1_down)
    print "number of entries that are up in WT2 is (test)", len (temp)
    assert len(wt1_down)==len(temp);
#QC
    se=cybert_collection()
    KOs=se.loaddb("./cybert_KOvsWT2cybert.db")
    #kos_up=se.find_differential_entries(se.samples[1])
    kos_up=se.find_differential_entries(se.samples[1],p=p);
    print "number of entries that are up in KO is", len (kos_up)
    kos_down=se.find_differential_entries(se.samples[0],p=p)
    print "number of entries that are DOWN in KO is", len (kos_down)
    sel=cybert_collection();
    OEs=sel.loaddb("./cybert_OEvsWT1cybert.db")
    oes_up=sel.find_differential_entries(sel.samples[1],p=p)
    print "number of entries that are up in OE is", len (oes_up)
    oes_down=sel.find_differential_entries(sel.samples[0],p=p)
    print "number of entries that are DOWN in OE is", len (oes_down)
    #out1=oes_up&kos_down;
    out1=kos_up-oes_up;
    #print "number of entries that are up in OE and down in KO is", len (out1)
    out2=oes_down-kos_down;

    pp=["lepr","sirt1"];
    m=set(sel.samples[0].entry_by_genesym2(pp));
    print "number of added additional entries are ", len(m);
    #out1=m;
    #out2=m;
    #out2=oes_up&kos_down-wt1_up;
    print "number of entries that are up in KO and not up in OE is", len (out1)
    print "number of entries that are DOWN in OE and not down in KO is", len (out2)
    #print "number of entries that are up in OE and down in KO but not up in WT1 compared to WT2 is", len (out2)
    #print "number of entries that are up in OE and up in KO is", len (oes_up&kos_up)
    #print "number of entries that are up in OE and up in KO but not up in WT1 compared to WT2 is", len (oes_up&kos_up-wt1_up)
    out3=oes_up-kos_up;
    print "number of entries that are up in OE and not up in KO is", len (out3)
    out4=kos_down-oes_down;
    print "number of entries that are down in KO and not down in OE", len (out4)

    #out3=m;
    #out4=m;
    #pvals(out1,"./csv/OE+KO-.csv")
    #pvals(out2,"./csv/OE+KO-NOWT.csv")
    pvals(out3,"./csv/OE+.csv")
    pvals(out4,"./csv/KO-.csv")
    
    pvals(out1,"./csv/KO+.csv")
    pvals(out2,"./csv/OE-.csv")
    #pvals(out3,"./csv/OE-KO+.csv")
    #pvals(out4,"./csv/OE-KO+NOWT.csv")
    print "These have a total intersection of {0:d}(KO+OE-), {1:d}(KO-OE+) and a total union of length {2:d}".\
    format(len(out1&out2),len(out3&out4),len(out1|out2|out3|out4));
def phos():
    #now we proceed to do something unrelate:
    se=cybert_collection()
    se.loaddb("./cybert_KOvsWT2cybert.db")

    sel=cybert_collection();
    sel.loaddb("./cybert_OEvsWT1cybert.db")

    header="\t".join(["genesym","IlluminaProbe","WT1","OE","OE+","WT1 BSD","OE BSD","pOE","WT2","KO","KO-","WT2 BSD","KO BSD","pKO"]);
    pp=["lep","lepr","sirt1","ptp1b","pten","ptpre","tcptp","ptpn11","ptprd","ptpro","ptprr","ptprf","ptpn5","ppn14","ptpn3","ptpn6"]
    WT1=sel.samples[0].entry_by_genesym2(pp);
    WT1=sorted(WT1,key=cybert_entry.__hash__);
    OE=sel.samples[1].entry_by_genesym2(pp);
    OE=sorted(OE,key=cybert_entry.__hash__);
    WT2=se.samples[0].entry_by_genesym2(pp);
    WT2=sorted(WT2,key=cybert_entry.__hash__);
    KO=se.samples[1].entry_by_genesym2(pp);
    KO=sorted(KO,key=cybert_entry.__hash__);
    assert len(WT1)==len(WT2);
    h=open("./txts/phosphatases.txt","w")
    h.write(header);
    h.write("\n")
    for i in xrange(len(WT1)):
        assert(WT1[i]==OE[i] and WT2[i]==KO[i])
        temp="\t".join([\
            WT1[i].genesym, WT1[i].refid, str(WT1[i].mean), str(OE[i].mean),str(OE[i].mean-WT1[i].mean),
            str(WT1[i].bsd), str(OE[i].bsd),str(WT1[i].pval),
            str(WT2[i].mean), str(KO[i].mean),str(KO[i].mean-WT2[i].mean),
            str(WT2[i].bsd), str(KO[i].bsd),str(WT2[i].pval)
            ]);
        h.write(temp);
        h.write("\n");
    h.close();
    return;

def pvals(listin,fname):
    f=open(fname,"w")
    for i in listin:
        f.write("\t".join([i.refid.upper()]))
        #f.write(str(i.pval));
        f.write("\n")
    f.close();
    return


def main():
    #c=cybert_collection();
    #sirt1()
#    f=filtered_collection()
#    f.load("GSE3142_highly_expressed_by_uniref.txt") 
#    f.savedb();
    #f.loaddb();
#    print f.entries
    #f.find_all_pairwise_unique()
#    f.find_all_unique_genes()
    #sirt1();
    #cc1=cybert_collection();
    #cc2=cybert_collection();
    #cc3=cybert_collection();
    #cc2.load("./txts/cybert_KOvsWT2.txt")
    #cc1.load("./txts/cybert_OEvsWT1.txt")
    #cc3.load("./txts/cybert_WT1vsWT2.txt")
    #cc1.load_gene_syms();
    #cc2.load_gene_syms();
    #cc3.load_gene_syms();
    #cc.loaddb();
    #cc2.plotpvals();
    #cc.find_all_unique_refs()
    #cc.find_all_unique_overexpressed_refs([entry.data for entry in f.entries]);
    #merged_ids=[f.entries[0].data|f.entries[1].data,\
    #        f.entries[2].data|f.entries[3].data,\
    #        f.entries[4].data|f.entries[5].data,\
    #        f.entries[6].data|f.entries[7].data,\
    #        f.entries[8].data|f.entries[9].data];
    #cc.find_all_unique_differential_refs(merged_ids);
    #cc.find_all_unique_differential_refs(merged_ids,True);
    #cc.find_all_unique_refs(merged_ids);
    #cc1.savedb();
    #cc2.savedb();
    #cc3.savedb();
    #sirt1();
    phos();
    return;
#main();
if __name__=="__main__":
    main();
