#!/home/baldig/shared_libraries/centos64/pkgs/python/2.6.5/bin/python
#little script to split up the sirt1 data

f=open("rawdata-processed.txt","r");
line = f.readline();
labels=[x.strip() for x in line.strip().split() if x.strip()!=""];
print labels;
indexes={};
for label in labels:
    indexes[label]=labels.index(label);
print indexes
KO_WT=["POOLA","POOLB","POOLC"];
KO_KO=["POOLD","POOLD","POOLE"];
#OE_WT=["A.AVG_Signal","B.AVG_Signal","C.AVG_Signal","D.AVG_Signal"];
OE_OE="ABCD";
OE_WT="EFGH";
#KO_WT+=OE_WT.split();
#print KO_WT;
#OE_WT=KO_WT;
#OE_WT="EFGH";
#real index lookups

#now we want to generate three pairwise conditioned data for cyber T, tab demilimited
g3=open("WT1vsWT2.txt","w");
g2=open("OEvsWT1.txt","w");
g1=open("KOvsWT2.txt","w");
g4=open("rownames.txt","w")
template="{0}\t{1}\t{2}\n"
temp=labels;
ko_wt=[temp[indexes[i]] for i in KO_WT];
ko_ko=[temp[indexes[i]] for i in KO_KO];
oe_oe=[temp[indexes[i]] for i in OE_OE];
oe_wt=[temp[indexes[i]] for i in OE_WT];
g1.write("#"+template.format("Probe",\
        "\t".join(["WT"+str(i+1) for i in xrange(len(ko_wt))]),\
        "\t".join(["KO"+str(i+1) for i in xrange(len(ko_ko))])))
g2.write("#"+template.format("Probe",\
        "\t".join(["WT"+str(i+1) for i in xrange(len(oe_wt))]),\
        "\t".join(["OE"+str(i+1) for i in xrange(len(oe_oe))])))
g3.write("#"+template.format("Probe",\
        "\t".join(["WT_KO"+str(i+1) for i in xrange(len(ko_wt))]),\
        "\t".join(["WT_OE"+str(i+1) for i in xrange(len(oe_wt))])))
for line in f:
    temp=line.strip().split("\t");
    assert(len(temp)==len(labels))
#note here some items are empty strings
    ko_wt=[temp[indexes[i]] for i in KO_WT];
    ko_ko=[temp[indexes[i]] for i in KO_KO];
    oe_oe=[temp[indexes[i]] for i in OE_OE];
    oe_wt=[temp[indexes[i]] for i in OE_WT];
    name=temp[indexes["ProbeID"]];
    g4.write(name+"\t"+temp[indexes["TargetID"]]);
    
    g1.write(template.format(name,\
            "\t".join([i for i in ko_wt]),\
            "\t".join([i for i in ko_ko])))
    g2.write(template.format(name,\
            "\t".join([i for i in oe_wt]),\
            "\t".join([i for i in oe_oe])))
    g3.write(template.format(name,\
            "\t".join([i for i in ko_wt]),\
            "\t".join([i for i in oe_wt])))
#done





