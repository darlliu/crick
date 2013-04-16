#specialized scripts for CD271 analysis
import math
CUTOFF=5e-3
col_template="data.addColumn('{0}','{1}');\n"
row_template="['{0}', new Date ({1:d},1,1), {2:f}, {3:d}, {4:f}, {5:f}, {6:f},'{7}'],\n"
def do_draw(targets,name):
    """
    For each target in the list draw the ones where pval is smaller than
    threashold and draw them depending on their diffmain
    """
    accum=set([]);
    for target in targets:
        #first accumulate the entries to plot
        UP=set([
            i for i in xrange(len(target.samples[0].entries))\
                    if target.samples[0].entries[i].d(target.samples[1].entries[i],CUTOFF)==-1
            ])
        DOWN=set([
            i for i in xrange(len(target.samples[0].entries))\
                    if target.samples[0].entries[i].d(target.samples[1].entries[i],CUTOFF)==1
            ])
        #relatively unsafe but fast
        print "For {0} we have {1:d} ups and {2:d} downs with cutoff pval of {3:.5f}".format(target.name,len(UP),len(DOWN),CUTOFF)
        accum|=UP;
        accum|=DOWN;
    print "total number of genes to plot {0:d}".format(len(accum))
    #now filter based on accum the entries, getting these info:
    #genesym, probeid, refseq (will be the order), mean, pval, (1/pval)*toone(diffmean)
    cols="";
    cols+=col_template.format("string","GeneSym")
    cols+=col_template.format("date","Date")
    cols+=col_template.format("number","Change")
    cols+=col_template.format("number","Order")
    cols+=col_template.format("number","mean_wt")
    cols+=col_template.format("number","mean_exp")
    cols+=col_template.format("number","pval")
    cols+=col_template.format("string","ProbeID")
    template=open("./motion_template.html","r").read()
    f=open("motionchart_"+name+".html","w")
    rows=""
    index=1900
    for target in targets:
        index+=1;
        WT=[target.samples[0].entries[i] for i in accum]
        KO=[target.samples[1].entries[i] for i in accum]
        WT=sorted(WT)
        KO=sorted(KO)
        assert len(WT)==len(KO)
        assert WT[0]==KO[0]
        assert WT[-1]==KO[-1]
        #we do some QC so not to mess up the order
        for wt,ko in zip(WT,KO):
            assert wt==ko;
            line=row_template.format(
                    wt.genesym, index, -math.log(1e-5+wt.pval) if float(wt.mean)>float(ko.mean) else math.log(1e-5+wt.pval),
                    WT.index(wt),wt.mean,ko.mean,wt.pval,wt.refid
                    )
            rows+=line;
    f.write(template.format(cols,rows))
    f.close();


def main(c):
    for col in c.cols:
        print col.name
    diffs=[
            "CD271_all_together0",
            #"CD271_one_one0","CD271_one_one1","CD271_one_one2","CD271_one_one3",
            "CD271_two_two_10","CD271_two_two_11",
            "CD271_two_two_20","CD271_two_two_21",
            "CD271_two_two_30","CD271_two_two_31"
            ]
    xenos=[
            "CD271_all_together1",
            #"CD271_one_one4",
            "CD271_two_two_12",
            "CD271_two_two_22",
            "CD271_two_two_32"
            ]
    older=[
            "CD271_all_together2",
            #"CD271_one_one4",
            "CD271_two_two_12",
            "CD271_two_two_22",
            "CD271_two_two_32"
            ]
    do_draw(filter(lambda x:x.name in diffs, c.cols),"Patients" )
    #do_draw(filter(lambda x:x.name in xenos, c.cols),"Xeno")
    #do_draw(filter(lambda x:x.name in older, c.cols), "Older")
    return;
if __name__=="__main__":
    main()
