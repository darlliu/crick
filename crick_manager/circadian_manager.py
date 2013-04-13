import os
import math
import cPickle as pk
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import heatmap as hp
class analyzed_entry (object):
    def __init__(self):
        self.uid=-1;
        self.disrupted=0;
        self.p=[];
        self.mdiff=[];
        self.mean=[];
        #contains both wt and ko so lol
        self.bsd=[];
        #same
        self.pjtk=[];
        self.cir=[];
        self.disrupts=[];
        self.lags=[];
        self.genesym="";
        self.geneid="";
        self.probeid="";
        self.timepts=[];
        self.flag=0;
    def check_disrupted(self, index=0,p=1e-7):
        if (self.p[index]>p):
            return 0;
        elif (self.mdiff[index]>0):
            return 1;
        #indicate mutant is overexpressed compared to WT
        else:
            return -1;
    def qc(self):
        if len(self.mdiff)==len(self.p):
            if len(self.mean)==len(self.bsd):
                if len(self.mean[0])==len(self.bsd[0]):
                    return;
        else:
            print len(self.mdiff),len(self.p),len(self.disrupts),len(self.cir), len(self.mean[0]);
            raise IndexError;
    def get_all_disruptions(self):
        return [str(self.check_disrupted(i)) for i in xrange(len(self.mdiff))];
    #def return_smart_sort_value(self):
        #"""give a value for sorting depending on the expression levels"""
        #def normalize(means):
            #mean=sum(means)/len(means);
            #means=[i/mean-1 for i in means];
            #return means[0];
        #out=0;
        #if self.flag==0:
            #out=self.mdiff;
            #out=normalize(out);
        #elif self.flag==1:
            #out=self.mean[0];
            #out=normalize(out);
        #elif self.flag==2:
            #out=self.mean[1];
            #out=normalize(out);
        #return out;
    #def __cmp__(self,another):
        #"""compares two entries"""
        #this=self.return_smart_sort_value();
        #that=another.return_smart_sort_value();
        #if this > that:
            #return 1;
        #elif this == that:
            #return 0;
        #else:
            #return -1;
    def check_WT_KO_circadian(self,index=0):
        identifier=0;
        if 'WT' in self.cir[index] and max(self.mean[0])/min(self.mean[0])>=1.5:
            identifier=1;
        elif 'KO' in self.cir[index] and max(self.mean[1])/min(self.mean[1])>=1.5:
            identifier=-1;
        elif 'BOTH' in self.cir[index] \
                and max(self.mean[0])/min(self.mean[0])>=1.5\
                and max(self.mean[1])/min(self.mean[1])>=1.5\
                :
            identifier=2;
        return identifier;
    def export_circadian_plot(self,prefix="default_",index=0, db=None, bypass=False):
        fig=plt.figure();
        wtmean=np.array(self.mean[0]);
        komean=np.array(self.mean[1]);
        wterr=np.array(self.bsd[0]);
        koerr=np.array(self.bsd[1]);
        print wtmean, wterr, komean, koerr;
        plt.errorbar(self.timepts,wtmean,xerr=0,yerr=wterr,label="WT");
        plt.errorbar(self.timepts,komean,xerr=0,yerr=koerr,label="KO");
        plt.title("Circadian behavior of gene: "+self.genesym);
        plt.xlim([self.timepts[0]-8, self.timepts[-1]+8]);
        plt.xlabel("Time, hours");
        plt.ylabel("Expression lv");
        plt.legend();
        ax=fig.add_subplot(111);
        if self.check_WT_KO_circadian(index)==2 or bypass:
            print "Both cir ",str( self.pjtk[0])[:len('%.*f' % (5, self.pjtk[0]))] , str( self.pjtk[-1])[:len('%.*f' % (5, self.pjtk[-1]))]
            pairxy=(self.timepts[0],wtmean[0]);
            ax.annotate('wt circadian p={0:.3e}'.format( self.pjtk[0]), xy=pairxy,xycoords='data',\
                xytext=(-50,30),textcoords='offset points',\
                arrowprops=dict(arrowstyle="->"));
            pairxy=(self.timepts[0],komean[0]);
            ax.annotate('ko circadian p={0:.3e}'.format( self.pjtk[-1]), xy=pairxy,xycoords='data',\
                xytext=(50,-30),textcoords='offset points',\
                arrowprops=dict(arrowstyle="->"));
        elif self.check_WT_KO_circadian(index)==1:
            print "WT circ", str( self.pjtk[0])[:len('%.*f' % (5, self.pjtk[0]))]
            pairxy=(self.timepts[0],wtmean[0]);
            ax.annotate('wt circadian p={0:.3e}'.format( self.pjtk[0]), xy=pairxy,xycoords='data',\
                xytext=(-50,30),textcoords='offset points',\
                arrowprops=dict(arrowstyle="->"));
        elif self.check_WT_KO_circadian(index)==-1:
            print "KO circ", str( self.pjtk[-1])[:len('%.*f' % (5, self.pjtk[-1]))]
            pairxy=(self.timepts[0],komean[0]);
            ax.annotate('ko circadian p={0:.3e}'.format( self.pjtk[-1]), xy=pairxy,xycoords='data',\
                xytext=(-50,30),textcoords='offset points',\
                arrowprops=dict(arrowstyle="->"));
        lookup=[".12",".18",".24",".30"];
        if db!=None:
            for i in xrange(len(self.timepts)):
                myy1=db["WT"+lookup[i]];
                myy2=db["KO"+lookup[i]];
                print myy1,wtmean, myy2, komean;
                myx=[self.timepts[i] for j in myy1];
                plt.plot(myx,myy1,marker='o')
                plt.plot(myx,myy2,marker='7')
        plt.savefig(prefix+str(self.uid)+".png");
        return

class crick_pipeline_result_analyzer(object):
    def __init__(self):
        self.entries=[];
        self.lookup={};
        self.timepoints=0;
        return;
    def load (self):
        try:
            f=open("Crick_table.xls","r");
        except IOError:
            print "Load error, try renaming input files!"
            raise;
        line=f.readline();
        temp=[x for x in line.strip().split() if x!=""];
        for j in range(0,len(temp)):
            self.lookup[temp[j]]=j;
        i=0;
        pvals=[x for x in temp if 'pVal' in x];
        amps=[x for x in temp if 'AMP' in x];
        #will be WT, KO amps
        times=[int (x.split('_')[-1]) for x in pvals]
        print times;
        diffmeans=[x for x in temp if 'diffMean_' in x];
        means1=[x for x in temp if 'MEAN_WT' in x];
        means2=[x for x in temp if 'MEAN_KO' in x];
        lags=[x for x in temp if '_LAG' in x];
        print lags
        bsd1=[x for x in temp if 'SD_WT' in x];
        bsd2=[x for x in temp if 'SD_KO' in x];
        circadian=[x for x in temp if 'circadian_' in x];
        disrupts=[x for x in temp if 'disrupted_' in x];
        pjtks=[x for x in temp if '_ADJP' in x];
        print pvals, diffmeans
        if len(pvals)!=len(diffmeans):
            print "Parsing length mismatch!!!";
            raise;
        else:
            self.timepoints=len(pvals);
        # get all pval and diffmean keys
        for line in f:
            if (i==0):
                i+=1;
                continue;
            temp=[x for x in line.strip().split() if x!=""];
            entry=analyzed_entry();
            entry.timepts=times;
            entry.uid=int(temp[self.lookup["uid"]]);
            entry.p=[float(temp[self.lookup[pkey]]) for pkey in pvals];
            entry.pjtk=[float(temp[self.lookup[pkey]]) for pkey in pjtks];
            entry.amps=[float(temp[self.lookup[pkey]]) for pkey in amps];
            entry.mdiff=[float(temp[self.lookup[key]]) for key in diffmeans];
            entry.mean.append([float(temp[self.lookup[key]]) for key in means1]);
            entry.mean.append([float(temp[self.lookup[key]]) for key in means2]);
            entry.bsd.append([float(temp[self.lookup[key]]) for key in bsd1]);
            entry.bsd.append([float(temp[self.lookup[key]]) for key in bsd2]);
            entry.lags=[float(temp[self.lookup[key]]) for key in lags];
            entry.cir=[temp[self.lookup[key]] for key in circadian];
            entry.disrupts=[int(temp[self.lookup[key]]) for key in disrupts];
            entry.geneid=temp[self.lookup["geneid"]];
            entry.probeid=temp[self.lookup["Probe_Set_ID"]];
            entry.genesym=temp[self.lookup["genesym"]];
            entry.qc();
            self.entries.append(entry);
        print len(self.entries), self.entries[0].p, self.entries[0].mdiff,self.entries[0].cir;
    def get_ups_and_downs(self, timept,p=1e-7):
        if timept>self.timepoints:
            raise;
        ups=[entry.uid for entry in self.entries if (entry.check_disrupted(timept,p)==1)];
        dns=[entry.uid for entry in self.entries if (entry.check_disrupted(timept,p)==-1)];
        return (ups,dns);
    def get_ups_and_downs2(self, timept,p=1e-7,foldchange=1.0):
        if timept>self.timepoints:
            raise;
        print "fold change is ", foldchange, "log is ", math.log(foldchange,2);
        ups=[entry.genesym.upper() for entry in self.entries\
                if ( (entry.mean[1][timept]-entry.mean[0][timept])>math.log(foldchange,2))];
        dns=[entry.genesym.upper() for entry in self.entries
                if ( (entry.mean[0][timept]-entry.mean[1][timept])>math.log(foldchange,2))];
        return (ups,dns);
    def get_ups_and_downs3(self, timept,p=1e-7,foldchange=1.0):
        if timept>self.timepoints:
            raise;
        print "fold change is ", foldchange, "log is ", math.log(foldchange,2);
        ups=[entry.genesym.upper() for entry in self.entries\
                if ( (entry.mean[1][timept]/entry.mean[0][timept])>foldchange)];
        dns=[entry.genesym.upper() for entry in self.entries
                if ( (entry.mean[0][timept]/entry.mean[1][timept])>foldchange)];
        return (ups,dns);
    def write_ups_and_downs(self,p=1e-7,foldchange=1.0, linear=False):
        print "fold change is ", foldchange
        fout=open("up_downs.txt","w");
        gout=open("overall_up_downs.txt",'w')
        for entry in self.entries:
            fout.write(str(entry.uid));
            fout.write('\t');
            fout.write("\t".join(entry.get_all_disruptions()));
            fout.write("\n");
        accum_ups=set([entry.uid for entry in self.entries]);
        accum_downs=set([entry.uid for entry in self.entries]);
        UP=[];
        DOWN=[];
        f=open("UP_GENESYMS.txt",'w')
        g=open("DOWN_GENESYMS.txt",'w')
        print "starting with", len(accum_ups),"entries";
        for i in  xrange (self.timepoints):
            tps=self.get_ups_and_downs(i,p);
            if (linear):
                tpg=self.get_ups_and_downs3(i,p, foldchange);
            else:
                tpg=self.get_ups_and_downs2(i,p, foldchange);

            f.write('\n'.join(tpg[0]))
            f.write('\n');
            g.write('\n'.join(tpg[1]))
            g.write('\n');
            gout.write(str(i)+"\t");
            print "up to time point:", i, "ups: ", len(tps[0]), "downs: ", len (tps[1]);
            gout.write(str(len(tps[0]))+"\t"+str(len(tps[1]))+"\t");
            UP.append(len(tpg[0]))
            DOWN.append(-len(tpg[1]))
            accum_ups&=set(tps[0]);
            accum_downs&=set(tps[1]);
            print "also, consecutive ups:", len (accum_ups), "consecutive downs: ", len(accum_downs);
            gout.write(str(len(accum_ups))+"\t"+str(len(accum_downs))+"\n");
        fout.close();
        gout.close();
        plt.figure();
        tp=[12,18,24,30];
        print UP, DOWN
        for i in xrange(self.timepoints):
            p1=plt.bar(tp[i],UP[i],color='r',label=str(UP[i]))
            p2=plt.bar(tp[i],DOWN[i],color='y',label=str(DOWN[i]))
            plt.legend((p1[0],p2[0]),('Ups','Downs'))
        plt.xticks(tp,[str(item) for item in tp])
        plt.text(20,min(DOWN)-min(DOWN)/10,"ups: "+str(UP))
        plt.text(20,min(DOWN)-min(DOWN)/5,"downs: "+str(DOWN))
        plt.text(20,min(DOWN)-min(DOWN)/2,"total "+str(len(self.entries)))
        plt.title("Number of up and down genes, fold>"+str(foldchange))
        plt.grid(which='major')
        plt.savefig(str(foldchange)+"fold_ups_and_downs.png")
        return;
    def write_circadians(self,index=0):
        #here we take p=0.05
        fout=open("circadians"+str(index)+".txt",'w');
        ref=open("lookup.cpkl.db","rb");
        DB=pk.load(ref)
        wts=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==1];
        kos=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==-1];
        both=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==2];
        print "WT circadian: ", len (wts), "KO circadians:", len(kos),\
                "both:",len(both),"neither:", len(self.entries)-len(wts)-len(kos)-len(both)
        fout.write("\t".join(['genesym','WT circadian','KO circadian','BOTH circadian','Disrupted time points']))
        fout.write("\t");
        fout.write("\t".join(['WT circadian P value', 'KO circadian P value']));
        fout.write("\t");
        fout.write("\t".join(['{0} mean levels at {1}'.format(i,j) for i in ['WT', 'KO']\
                for j in ['12 hours','18 hours','24 hours','30 hours']]))
        fout.write("\t");
        fout.write("\t".join(['{0} standard deviation (bayesian) at {1}'.format(i,j) for i in ['WT', 'KO']\
                for j in ['12 hours','18 hours','24 hours','30 hours']]))
        fout.write("\t");
        fout.write("\t".join(['p value for disrupted at {0}'.format(j)\
                for j in ['12 hours','18 hours','24 hours','30 hours']]))
        fout.write("\t");
        fout.write("\t".join(['WT Rep{0}, {1} hrs'.format(i,j)\
                for j in ['12 hours','18 hours','24 hours','30 hours']\
                for i in "123"\
                ]))
        fout.write("\t");
        fout.write("\t".join(['KO Rep{0}, {1} hrs'.format(i,j)\
                for j in ['12 hours','18 hours','24 hours','30 hours']\
                for i in "123"\
                ]))
        fout.write("\n");
        gout=open(str(index)+"circadians_ids_wt.txt",'w');
        for entry in wts:
            db=DB[entry.geneid];
            fout.write("\t".join([entry.genesym, 'Yes', 'No','No'\
                    ,str(entry.disrupts[index])]))
            fout.write("\t");
            fout.write("\t".join([str(item) for item in\
                    (entry.pjtk+entry.mean[0]+entry.mean[1]+\
                    entry.bsd[0]+entry.bsd[1]+entry.p\
                    +db["WT.12"]+db["WT.18"]+db["WT.24"]+db["WT.30"]\
                    +db["KO.12"]+db["KO.18"]+db["KO.24"]+db["KO.30"]\
                    )]));
            gout.write(entry.genesym);
            gout.write("\n");
            fout.write("\n");
        gout=open(str(index)+"circadians_ids_ko.txt",'w');
        for entry in kos:
            db=DB[entry.geneid];
            fout.write("\t".join([entry.genesym, 'No', 'Yes','No'\
                    ,str(entry.disrupts[index])]))
            fout.write("\t");
            fout.write("\t".join([str(item) for item in\
                    (entry.pjtk+entry.mean[0]+entry.mean[1]+\
                    entry.bsd[0]+entry.bsd[1]+entry.p+\
                    db["WT.12"]+db["WT.18"]+db["WT.24"]+db["WT.30"]\
                    +db["KO.12"]+db["KO.18"]+db["KO.24"]+db["KO.30"]\
                    )]));
            gout.write(entry.genesym);
            gout.write("\n");
            fout.write("\n");
        gout=open(str(index)+"circadians_ids_both.txt",'w');
        for entry in both:
            db=DB[entry.geneid];
            fout.write("\t".join([entry.genesym, 'Yes', 'Yes','Yes'\
                    ,str(entry.disrupts[index])]))
            fout.write("\t");
            fout.write("\t".join([str(item) for item in\
                    (entry.pjtk+entry.mean[0]+entry.mean[1]+\
                    entry.bsd[0]+entry.bsd[1]+entry.p+\
                    db["WT.12"]+db["WT.18"]+db["WT.24"]+db["WT.30"]\
                    +db["KO.12"]+db["KO.18"]+db["KO.24"]+db["KO.30"]\
                    )]));
            gout.write(entry.genesym);
            gout.write("\n");
            fout.write("\n");
        gout=open(str(index)+"circadians_all.txt",'w');
        for entry in both+wts+kos:
            gout.write(entry.genesym);
            gout.write("\n");
        fout.close();
        return
    def export_all_circadian_figs(self,index=0):
        subdir="./exports/img/";
        ref=open("lookup.cpkl.db","rb");
        DB=pk.load(ref)
        for entry in self.entries:
            flag=entry.check_WT_KO_circadian(index);
            if flag==0:
                continue;
            else:
                db=DB[entry.geneid]
                entry.export_circadian_plot(subdir,index,db);
    def export_all_circadian_figs_with_flags(self,index=0,flag=0):
        subdir="./exports/img/";
        converts={1:0,-1:1,2:2};
        flag=converts[flag];
        labels=["WT only", "KO only", "both"];
        for entry in self.entries:
            if entry.check_WT_KO_circadian()==0:
                continue;
            elif entry.check_WT_KO_circadian()==flag:
                entry.export_circadian_plot(subdir+labels(flag));
    

    def generate_table_elements(self,collection, linkout=False):
        out="";
        i=0;
        while i < len(collection):
            out+="<tr>\n";
            for j in xrange(5):
                if linkout:
                    out+="<td>\n<a href=\"./img/{0}.png\">\n".format(str(collection[i].uid));
                else:
                    out+="<td>\n<a href=\"#{0}\">\n".format(str(collection[i].uid));
                out+=collection[i].genesym;
                i+=1;
                out+="</td>\n</a >\n";
                if i>=len(collection):
                    break;
            out+="</tr>\n";
        return out;
    def generate_circadian_report(self,index=0,linkout=False):
        WTs=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==1];
        KOs=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==-1];
        BOTHs=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==2];
        ALLs=WTs+KOs+BOTHs;
        self.write_circadian_report(index,linkout,WTs,KOs,BOTHs,ALLs);
        return
    def generate_core_gene_report(self,index=3,linkout=False):
        coref=open("coregenes.list","r");
        temp=[];
        ref=open("lookup.cpkl.db","rb");
        DB=pk.load(ref)
        for line in coref:
            temp.append(line.strip());#cautious grab
        cores=[item.strip().lower() for item in temp if item.strip()!=""];
        print "here are the core genes:", cores;
        into=[entry for entry in self.entries if entry.genesym.lower() in cores];
        print len(into)
        #here we do reverse lookup to avoid missing
        for stuff in into:
            db1=DB[stuff.geneid];
            stuff.export_circadian_plot(prefix="./exports/img/",index=2, db=db1,bypass=True);
        self.write_circadian_report(index,linkout,into,[],[], into,templatename="template_core.html");
        return;

    def write_circadian_report(self,index=0,linkout=False,\
            WTs=[], KOs=[],BOTHs=[],ALLs=[],\
            templatename="template.html"):
        """generated a more prettified report based on info parsing,
        uses bootstrap from twitter."""
        imgdir="./img/";
        htmldir="./exports/";
        template=open(templatename,'r');
        #first, get a list of all circadian wts and kos and boths
        print "alls", len (ALLs);
        outstr="";
        tempstr="";
        block="";
        images="";
        spacer="""<div class="span8 offset2"></div>""";
        table1=self.generate_table_elements(WTs,linkout);
        table2=self.generate_table_elements(KOs,linkout);
        table3=self.generate_table_elements(BOTHs,linkout);
        print "tables", len (table1);
        #grab block template for images
        flag=False;
        for line in template:
            if "%%block%%" in line:
                flag=~flag;
            while flag:
                block+=line;
                break;
        #write out all images inserts
        template.close();
        print "blocks", len (block);
        i=0;
        while i!=len(ALLs):
            temp=block.replace("$NAME1$",str(ALLs[i].uid));
            temp=temp.replace("$SRC1$",imgdir+str(ALLs[i].uid)+".png");
            i+=1;
            if i>=len(ALLs):
                break;
            temp=temp.replace("$NAME2$",str(ALLs[i].uid));
            temp=temp.replace("$SRC2$",imgdir+str(ALLs[i].uid)+".png");
            i+=1;
            images+=temp;
        if linkout:
            images="";
        print "images", len (images);
        flag1=flag2=flag3=flag=False;
        template=open(templatename,'r');
        for line in template:
            if "%%table1%%" in line:
                flag1=~flag1;
            elif "%%table2%%" in line:
                flag2=~flag2;
            elif "%%table3%%" in line:
                flag3=~flag3;
            elif "%%block%%" in line:
               flag=~flag;
            while flag1:
                outstr+=table1;
                table1=spacer;
                break;
            while flag2:
                outstr+=table2;
                table2=spacer;
                break;
            while flag3:
                outstr+=table3;
                table3=spacer;
                break;
            while flag:
                outstr+=images;
                images=spacer;
                break;
            if ~(flag and flag1 and flag2 and flag3):
                outstr+=line;
        labels=['0.001','0.01','0.05','coregene']
        if linkout:
            fout=open("./exports/report_linkout_{0}.html".format(labels[index]),"w");
        else:
            fout=open("./exports/report_big_{0}.html".format(labels[index]),"w");
        fout.write(outstr);
    
    def generate_pie_chart(self):
        """generate a simple pie chart that describes the circadian behavior"""
        plt.figure();
        pvals=[0.001,0.01,0.05];
        for index in [0, 1 ,2]:
            wts=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==1];
            kos=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==-1];
            both=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==2];
            X=[len(wts),len(kos),len(both)];
            plt.subplot(311+index);
            plt.pie(X,labels=["WT: "+str(len(wts)),"KO: "+str(len(kos)),"BOTH: "+str(len(both))]);
            plt.title("number of circadian genes, cutoff p="+str(pvals[index]));
        plt.savefig("piecharts.png");    
        return;
    def generate_heat_map_one_by_one(self,index=0,flag=True,tosort=1):
        """generate a heat map of circadian and non-circadian genes expression profiles"""
        if flag:
            circs=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)!=0];
        else:
            circs=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==0];
        for entry in circs:
            entry.flag=tosort;
        circs=sorted(circs);
        x=[12,18,24,30];
        y=xrange(len(circs));
        X,Y=np.meshgrid(x,y);
        labels={0:"df_sorted_", 1:"WT_sorted_",2:"KO_sorted_"};
        plt.figure();
        plt.title("heat map of expression levels, side by side (WT,KO,WT,KO...)");
        plt.subplot(121);
        plt.title("WT and KO, side by side");
        #for entries in circs:
        #    entries.flag=1;
        #mean_wt=[entry.mean[0] for entry in circs];
        mean_wt=[];
        for entry in circs:
            temp=[];
            for i in xrange(4):
                temp.append(entry.mean[0][i]);
                temp.append(entry.mean[1][i]);
            mean_wt.append(temp)
        Z=np.array([[mean_wt[j][i] for i in xrange(8)] for j in y]);
        plt.pcolormesh(Z,vmax=14, vmin=0);
        plt.colorbar();
        plt.grid(True, which='major');
        plt.xticks([0,2,4,6],['12','18','24','30'])
        #frame1 = plt.gca();
        #frame1.axes.get_xaxis().set_ticks([0,2,4,6],['12','18','24','30'])
        plt.xlabel("hours");
        plt.ylabel("genes");


        plt.subplot(122);
        plt.title("KO-WT");
        mean_df=[entry.mdiff for entry in circs];
        Z=np.array([[mean_df[j][i] for i in xrange(4)] for j in y]);
        plt.pcolormesh(Z,vmax=1,vmin=-1);
        plt.colorbar();
        plt.grid(True, which='major');
        plt.xticks([0,1,2,3],['12','18','24','30'])
        plt.xlabel("hours");
        #toplot=[[float(i),float(j),meanwt[k]] for k in xrange(4) for i in [12,18,24,30] for j in xrange(len(circs)) for meanwt in mean_wt];
        #pp=np.array(toplot);
        #plt.imshow(toplot,cmap=cm.RdYlGn, interpolation='nearest', extent=[12,30, 1, len(circs)], \
        #        vmax=15, vmin=5);
        #plt.hexbin(np.array(X),np.array(Y),Z);
        if flag:
            plt.savefig(labels[tosort]+"circadian_heat.png");
        else:
            plt.savefig(labels[tosort]+"non_circadian_heat.png");
    def generate_heat_map(self,index=0,flag=True,tosort=1):
        """generate a heat map of circadian and non-circadian genes expression profiles"""
        def normalize(means):
            mm=max([abs(i) for i in means]);
            mi=min([abs(i) for i in means]);
            means=[(i-mi)/(mm-mi)-0.5 for i in means];
            return means;
        if flag:
            circs=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)!=0];
        else:
            circs=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==0];
        for entry in circs:
            entry.flag=tosort;
        circs=sorted(circs,key=lambda a:-a.lags[tosort]);
        x=[12,18,24,30,12,18,24,30];
        y=xrange(len(circs));
        X,Y=np.meshgrid(x,y);
        labels={-1:"df_sorted_", 0:"WT_sorted_",1:"KO_sorted_"};
        plt.figure();
        plt.title("heat map of expression levels");
        plt.subplot(121);
        plt.title("WT and then KO");
        #for entries in circs:
        #    entries.flag=1;
        means=[normalize(entry.mean[0])+normalize(entry.mean[1]) for entry in circs];
        Z=np.array([[means[j][i] for i in xrange(8)] for j in y]);
        plt.pcolormesh(Z,vmax=0.5, vmin=-0.5);
        plt.colorbar();
        plt.grid(True, which='major');
        plt.xticks([0,1,2,3,4,5,6,7],['12','18','24','30','12','18','24','30'])
        plt.xlabel("hours");
        plt.ylabel("genes");


        plt.subplot(122);
        plt.title("KO-WT");
        mean_df=[normalize(entry.mdiff) for entry in circs];
        Z=np.array([[mean_df[j][i] for i in xrange(4)] for j in y]);
        plt.pcolormesh(Z,vmax=1,vmin=-1);
        plt.grid(True, which='major');
        plt.xticks([0,1,2,3],['12','18','24','30'])
        frame1 = plt.gca();
        frame1.axes.get_yaxis().set_ticks([])
        plt.colorbar();
        plt.xlabel("hours");
        #toplot=[[float(i),float(j),meanwt[k]] for k in xrange(4) for i in [12,18,24,30] for j in xrange(len(circs)) for meanwt in mean_wt];
        #pp=np.array(toplot);
        #plt.imshow(toplot,cmap=cm.RdYlGn, interpolation='nearest', extent=[12,30, 1, len(circs)], \
        #        vmax=15, vmin=5);
        #plt.hexbin(np.array(X),np.array(Y),Z);
        pp=['0.001','0.01','0.05'];
        if flag:
            plt.savefig(labels[tosort]+"circadian_heat"+"p_"+str(pp[index])+".png");
        else:
            plt.savefig(labels[tosort]+"non_circadian_heat"+"p_"+str(pp[index])+".png");
    def generate_heat_map3(self,index=0,flag=True,tosort=1):
        """generate a heat map of circadian and non-circadian genes expression profiles"""
        if flag:
            circs=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)!=0];
        else:
            circs=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==0];
        for entry in circs:
            entry.flag=tosort;
        circs=sorted(circs);
        x=[12,18,24,30];
        y=xrange(len(circs));
        X,Y=np.meshgrid(x,y);
        labels={0:"df_sorted_", 1:"WT_sorted_",2:"KO_sorted_"};
        plt.figure();
        plt.title("heat map of expression levels");
        plt.subplot(131);
        plt.title("WT");
        #for entries in circs:
        #    entries.flag=1;
        mean_wt=[entry.mean[0] for entry in circs];
        Z=np.array([[mean_wt[j][i] for i in xrange(4)] for j in y]);
        plt.pcolormesh(Z,vmax=14, vmin=0);
        plt.colorbar();
        plt.grid(True, which='major');
        plt.xticks([0,1,2,3],['12','18','24','30'])
        plt.xlabel("hours");
        plt.ylabel("genes");

        plt.subplot(132);
        plt.title("KO");
        #for entry in circs:
        #    entry.flag=2;
        #circs=sorted(circs);
        mean_ko=[entry.mean[1] for entry in circs];
        Z=np.array([[mean_ko[j][i] for i in xrange(4)] for j in y]);
        plt.pcolormesh(Z,vmax=14,vmin=0);
        plt.grid(True, which='major');
        plt.xticks([0,1,2,3],['12','18','24','30'])
        plt.colorbar();
        frame1 = plt.gca();
        frame1.axes.get_yaxis().set_ticks([])
        plt.xlabel("hours");

        plt.subplot(133);
        plt.title("KO-WT");
        mean_df=[entry.mdiff for entry in circs];
        Z=np.array([[mean_df[j][i] for i in xrange(4)] for j in y]);
        plt.pcolormesh(Z,vmax=1,vmin=-1);
        plt.grid(True, which='major');
        plt.xticks([0,1,2,3],['12','18','24','30'])
        frame1 = plt.gca();
        frame1.axes.get_yaxis().set_ticks([])
        plt.colorbar();
        plt.xlabel("hours");
        #toplot=[[float(i),float(j),meanwt[k]] for k in xrange(4) for i in [12,18,24,30] for j in xrange(len(circs)) for meanwt in mean_wt];
        #pp=np.array(toplot);
        #plt.imshow(toplot,cmap=cm.RdYlGn, interpolation='nearest', extent=[12,30, 1, len(circs)], \
        #        vmax=15, vmin=5);
        #plt.hexbin(np.array(X),np.array(Y),Z);
        pp=['0.001','0.01','0.05'];
        if flag:
            plt.savefig(labels[tosort]+"circadian_heat"+"p_"+str(pp[index])+".png");
        else:
            plt.savefig(labels[tosort]+"non_circadian_heat"+"p_"+str(pp[index])+".png");
    def generate_heat_new(self,index=1):
        circs=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)!=0];
        noncircs=[entry for entry in self.entries if entry.check_WT_KO_circadian(index)==0];
        mean_wt=[entry.mean[0] for entry in circs];
        mean_ko=[entry.mean[1] for entry in circs];
        hp.heatmap(np.array(mean_wt),"WT");
        hp.heatmap(np.array(mean_ko),"KO");
        return
def main():
    c=crick_pipeline_result_analyzer();
    c.load();
    c.write_ups_and_downs();
    #c.write_circadians(0);
    c.write_circadians(1);
    #c.write_circadians(2);
    index=1;
    c.export_all_circadian_figs(index);
    #c.generate_circadian_report(index);
    c.generate_circadian_report(index,linkout=True);
    #index=1;
    #c.generate_circadian_report(index);
    #c.generate_circadian_report(index,linkout=True);
    #index=0;
    #c.generate_circadian_report(index);
    #c.generate_circadian_report(index,linkout=True);
    #c.generate_heat_map(index=2,flag=True,tosort=1);
    #c.generate_heat_map(index=2,flag=False,tosort=1);
    #c.generate_heat_map(index=2,flag=True,tosort=0);
    #c.generate_heat_map(index=2,flag=False,tosort=0);

    c.generate_heat_map(index=1,flag=True,tosort=1);
    c.generate_heat_map(index=1,flag=False,tosort=1);
    c.generate_heat_map(index=1,flag=True,tosort=0);
    c.generate_heat_map(index=1,flag=False,tosort=0);

    #c.generate_heat_map(index=0,flag=True,tosort=1);
    #c.generate_heat_map(index=0,flag=False,tosort=1);
    #c.generate_heat_map(index=0,flag=True,tosort=0);
    #c.generate_heat_map(index=0,flag=False,tosort=0);

    #c.write_ups_and_downs(foldchange=1.2,linear=True);
    #c.write_ups_and_downs(foldchange=1.3,linear=True);
    #c.write_ups_and_downs(foldchange=1.4,linear=True);
    #c.write_ups_and_downs(foldchange=1.5,linear=True);
    #c.write_ups_and_downs(foldchange=1.6,linear=True);
    #c.write_ups_and_downs(foldchange=1.7,linear=True);
    #c.write_ups_and_downs(foldchange=1.8,linear=True);
    #c.write_ups_and_downs(foldchange=1.9,linear=True);
    #c.write_ups_and_downs(foldchange=2.0,linear=True);
    #c.generate_heat_map_one_by_one(index=1,flag=True,tosort=1);
    #c.generate_heat_map_one_by_one(index=1,flag=False,tosort=1);
    #c.generate_heat_map_one_by_one(index=1,flag=True,tosort=0);
    #c.generate_heat_map_one_by_one(index=1,flag=False,tosort=0);
    #c.generate_heat_map_one_by_one(index=1,flag=True,tosort=2);
    #c.generate_heat_map_one_by_one(index=1,flag=False,tosort=2);

    c.generate_core_gene_report();
    #c.generate_pie_chart();
main();
               

        

