#specialized scripts for CD271 analysis
import math
import string
import os;
col_template="data.addColumn('{0}','{1}');\n"
row_template="['{0}', new Date ({1:d},1,1), {2:f}, {3:d}, {4:f}, {5:f}, {6:f},'{7}'],\n"
def do_draw(targets,name):
    """
    For each target in the list draw the ones where pval is smaller than
    threashold and draw them depending on their diffmain
    """
    CUTOFF=5e-3
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
                    wt.genesym, index, math.log(1e-5+wt.pval) if float(wt.mean)>float(ko.mean) else -math.log(1e-5+wt.pval),
                    WT.index(wt),wt.mean,ko.mean,wt.pval,wt.refid
                    )
            rows+=line;
    f.write(template.format(cols,rows))
    f.close();

def pvals(targets):
    """P value analysis on all pairwise data. Select genes whose """
    global unique_up, unique_down, accum_up,accum_down, five_ups, five_downs, UP1, DOWN1
    CUTOFF=5e-2
    curdir=os.getcwd()
    try:
        os.mkdir(str(CUTOFF))
    except:
        pass;
    os.chdir(str(CUTOFF))
    unique_up=set([i for i in xrange (len(targets[0].samples[0].entries))] );
    unique_down=set([i for i in xrange (len(targets[0].samples[0].entries))] );
    accum_up=set([]);
    accum_down=set([]);
    #find by indices -- fast

    target=targets[0];
    UP=set([
        target.samples[0].entries[i].genesym for i in xrange(len(target.samples[0].entries))\
                if target.samples[0].entries[i].d(target.samples[1].entries[i],CUTOFF)==-1
        ])
    DOWN=set([
        target.samples[0].entries[i].genesym for i in xrange(len(target.samples[0].entries))\
                if target.samples[0].entries[i].d(target.samples[1].entries[i],CUTOFF)==1
        ])
    UP1=set([
        i for i in xrange(len(target.samples[0].entries))\
                if target.samples[0].entries[i].d(target.samples[1].entries[i],CUTOFF)==-1
        ])
    DOWN1=set([
        i for i in xrange(len(target.samples[0].entries))\
                if target.samples[0].entries[i].d(target.samples[1].entries[i],CUTOFF)==1
        ])
    open("Allcombined_ups.txt","w").write("\n".join(UP))
    open("Allcombined_downs.txt","w").write("\n".join(DOWN))
    logs=open("stats.txt","w")
    for target in targets[1:]:
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
        unique_up&=UP;
        unique_down&=DOWN;
        accum_up|=UP;
        accum_down|=DOWN;
        logs.write( "Af target {0}  we have {1:d} ups and {2:d} \
                downs with cutoff pval of {3:.5f}".format(target.name,len(UP),len(DOWN),CUTOFF)
                )
        logs.write("\n")
        logs.write( "After intersecting {0}  we have {1:d} unique ups and {2:d} \
                unique downs with cutoff pval of {3:.5f}".format(target.name,len(unique_up),len(unique_down),CUTOFF)
                )
        logs.write("\n")
        logs.write("As for union  we have {1:d} ups and {2:d} \
                downs with cutoff pval of {3:.5f}".format(target.name,len(accum_up),len(accum_down),CUTOFF)
                )
        logs.write("\n")
        logs.write("\n")
        logs.write("\n")
    f1=open("unique_ups.txt","w");
    f2=open("unique_downs.txt","w");
    for i in unique_up:
        f1.write(target.samples[0].entries[i].genesym)
        f1.write("\n")
    for i in unique_down:
        f2.write(target.samples[0].entries[i].genesym)
        f2.write("\n")
    f1.close()
    f2.close()
    g1=open("five_ups.txt","w")
    g2=open("five_downs.txt","w")
    five_ups=set([]);
    five_downs=set([]);

    for i in accum_up:
        sym=targets[0].samples[0].entries[i].genesym
        temp=0;
        for j in targets:
            if j.samples[0].entries[i].d(j.samples[1].entries[i],CUTOFF)==-1:
                temp+=1;
        assert(temp>=1)
        print "up for {0} there are {1:d} over threshold targets"\
                .format(sym,temp)
        if temp>=5:
            five_ups|=set([i])
            g1.write(sym)
            g1.write("\n")

    for i in accum_down:
        sym=targets[0].samples[0].entries[i].genesym
        temp=0;
        for j in targets:
            if j.samples[0].entries[i].d(j.samples[1].entries[i],CUTOFF)==1:
                temp+=1;
        assert(temp>=1)
        print "down: for {0} there are {1:d} over threshold targets"\
                .format(sym,temp)
        if temp>=5:
            five_downs|=set([i])
            g2.write(sym)
            g2.write("\n")
    g1.close()
    g2.close()
    os.chdir(curdir)
    return



def folds(targets):
    """Fold change analysis on all single patients"""

def export(ones, couples, indices, name):
    """do export routine given raw data ones couples, indices to iterate over and output name"""
    curdir=os.getcwd()
    templatedir=curdir+'/template'
    tabletemp=string.Template(open(templatedir+"/table_template.html","r").read())
    linetemp=string.Template(open(templatedir+"/lineplot_template.html","r").read())
    indextemp=string.Template(open(templatedir+"/index_template.html","r").read())
    heattemp=string.Template(open(templatedir+"/heatmap_template.html","r").read())
    try:
        os.mkdir(name);
    except:
        pass
    os.chdir(name)
    try:
        os.mkdir("lines");
    except:
        pass
    #lines
    cols=rows=""
    col="data.addColumn('{0}','{1}');\n"
    cols+=col.format("number","N- expression")
    cols+=col.format("string","expname")
    cols+=col.format("string","expdetail")
    cols+=col.format("number","N+ expression")
    cols+=col.format("string","expname2")
    cols+=col.format("string","expdetail2")
    #these are presumed undefined
    row="[new Date(2010,1,{0:d}), {1:f},'{3}','{4}',{2:f},undefined, undefined],\n"
    for i in indices:
        inc=0
        rows=""
        for j in ones+couples:
            inc+=1;
            rows+=row.format(
                    inc, j.samples[0].entries[i].mean,j.samples[1].entries[i].mean,
                    j.name, "experiment/analysis"
                    )
        f=open("lines/"+j.samples[0].entries[i].genesym+"_lines.html","w")
        f.write(
                linetemp.substitute({"cols":cols,"rows":rows})
                )
        f.close()
    #generate tables
    cols=rows=""
    row="[{0}],"
    cols+=col.format("string","probe id")
    cols+=col.format("string","gene sym")
    for some in ones+couples:
        cols+=col.format("number","{0} N- mean".format(some.name))
        cols+=col.format("number","{0} N+ mean".format(some.name))
        cols+=col.format("number","{0} fold change".format(some.name))
        cols+=col.format("number","{0} p value".format(some.name))
        cols+=col.format("number","{0} BH".format(some.name))
    for i in indices:
        sym=ones[0].samples[0].entries[i].genesym
        probe=ones[0].samples[0].entries[i].refid
        nums=[];
        for some in ones+couples:
            nums+=[some.samples[0].entries[i].mean]
            nums+=[some.samples[1].entries[i].mean]
            nums+=[2**(nums[1]-nums[0])]
            nums+=[some.samples[0].entries[i].pval]
            nums+=[some.samples[0].entries[i].bh]
        out="'{0}','{1}',".format(sym,probe)
        out+=','.join([str(i) for i in nums]);
        rows+=row.format(out)
    open("tables.html","w").write(
            tabletemp.substitute({"cols":cols,"rows":rows})
            )




    #heatmaps 1
    cols=rows=""
    row0="data.setCell({0:d},{1:d},'{2}');\n"
    row="data.setCell({0:d},{1:d},{2:f});\n"
    cols+=col.format("string","genesym")
    for one in ones:
        cols+=col.format("number",one.name)
    inr=0;
    for i in indices:
        inc=0;
        rows+=row0.format(
                inr,inc, ones[0].samples[0].entries[i].genesym
                )
        for j in ones:
            inc+=1;
            #sign=1 if j.samples[0].entries[i].mean<j.samples[1].entries[i].mean else -1;
            rows+=row.format(
                inr,inc, (j.samples[1].entries[i].mean-j.samples[0].entries[i].mean)
                    )
        inr+=1;
    open("heatmap1.html","w").write(
            heattemp.substitute(
                {"cols":cols,"rows":rows,"which":"","rowcounts":str(len(indices))}
                ))



    cols=rows=""
    row0="data.setCell({0:d},{1:d},'{2}');\n"
    row="data.setCell({0:d},{1:d},{2:f});\n"
    cols+=col.format("string","genesym")
    for one in couples:
        cols+=col.format("number",one.name)
    inr=0;
    for i in indices:
        inc=0;
        rows+=row0.format(
                inr,inc, ones[0].samples[0].entries[i].genesym
                )
        for j in couples:
            inc+=1;
            sign=1 if j.samples[0].entries[i].mean<j.samples[1].entries[i].mean else -1;
            rows+=row.format(
                inr,inc, sign*math.log(j.samples[0].entries[i].pval)
                    )
        inr+=1;
    open("heatmap2.html","w").write(
            heattemp.substitute(
                {"cols":cols,"rows":rows,"rowcounts":str(len(indices)),"which":""}
                ))


    #generate menus
    col="<tr> {0}</td>\n"
    row="""<td> <a href="{0}"> {1} </a> </td>"""
    head="<th> {0} </th> "
    contents=rows="";
    cnt=0;
    for i in indices:
        sym=ones[0].samples[0].entries[i].genesym
        rows+=row.format(
                "lines/"+sym+"_lines.html",sym
                )
        cnt+=1;
        if cnt>10:
            cnt=0;
            contents+=col.format(rows)
            rows=""
    contents+=col.format(rows)
    open("index.html","w").write(
            indextemp.substitute(
                {"menu_content":contents, "menu_header": "", "title":name}
                )
            )

    os.chdir(curdir)
    return;



def main(c):
    for col in c.cols:
        print col.name
    folds=[
            "CD271_one_one0","CD271_one_one1","CD271_one_one2","CD271_one_one3"
            ]
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
    #do_draw(filter(lambda x:x.name in diffs, c.cols),"Patients" )
    pvals(filter(lambda x:x.name in diffs, c.cols) )
    export(filter(lambda x:x.name in folds, c.cols),filter(lambda x:x.name in diffs, c.cols) , unique_up, "Unique_Up")
    export(filter(lambda x:x.name in folds, c.cols),filter(lambda x:x.name in diffs, c.cols) , unique_down, "Unique_Down")
    export(filter(lambda x:x.name in folds, c.cols),filter(lambda x:x.name in diffs, c.cols) , five_ups, "five_ups")
    export(filter(lambda x:x.name in folds, c.cols),filter(lambda x:x.name in diffs, c.cols) , five_downs, "five_downs")
    export(filter(lambda x:x.name in folds, c.cols),filter(lambda x:x.name in diffs, c.cols) , UP1, "Allcombined_ups")
    export(filter(lambda x:x.name in folds, c.cols),filter(lambda x:x.name in diffs, c.cols) , DOWN1, "Allcombined_downs")
    #do_draw(filter(lambda x:x.name in xenos, c.cols),"Xeno")
    #do_draw(filter(lambda x:x.name in older, c.cols), "Older")
    return;
if __name__=="__main__":
    main()
