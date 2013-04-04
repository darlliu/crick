#!/home/baldig/shared_libraries/centos64/pkgs/python/2.6.5/bin/python
import cPickle as pickle

f=open("./MouseWG-6_V2_0_R3_11278593_A.txt")
h=open("./IlluminaAllRefs.txt",'w')
for line in f:
    header=[i.strip() for i in line.split("\t") if i.strip()!=""];
    break;
lookup={};
for i in header:
    lookup[i]=header.index(i)

print header, lookup
DB={}
for line in f:
    temp=[i.strip() for i in line.split("\t")];
    if len (temp)!=len(header):
        continue;
    #print temp[lookup["Search_Key"]]
    #print temp[lookup["Probe_Id"]]
    content={};
    for i in header:
        content[i]=temp[lookup[i]];
    DB[temp[lookup["Probe_Id"]]]=content;
    h.write(temp[lookup["Probe_Id"]]);
    h.write("\n");
print "Grabbed " , len(DB), "entries"
g=open("IlluminaProbeLookup.pkl","wb")
pickle.dump(DB,g);
print "done pickling"
g.close();
