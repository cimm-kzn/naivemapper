from CGRtools.CGRcore import CGRcore
from CGRtools.FEAR import FEAR
from CGRtools.RDFread import RDFread
from CGRtools.SDFwrite import SDFwrite

CGR=CGRcore(type='0', stereo=False, b_templates=None, balance=0, c_rules= None, e_rules=None)
fear=FEAR()
SDF=SDFwrite(open("out_fold2.sdf", 'w'))

with open("output_fold2pred.rdf") as predfile, open("output_fold2test.rdf") as testfile:
    ok=0
    nok=0
    for pred,test in zip(RDFread(predfile).readdata(),RDFread(testfile).readdata()):
        predCGR=CGR.getCGR(pred)
        testCGR=CGR.getCGR(test)
        predHash=fear.getreactionhash(predCGR)
        testHash=fear.getreactionhash(testCGR)
        #print("pred:" + predHash)
        #print("test:"+testHash)
        if predHash==testHash:
            ok+=1
            print("O'k")
        else:
            nok+=1
            print("Not O'k")
        tmp=CGR.getformattedcgr(predCGR)
        tmp['meta']={}
        SDF.writedata(tmp)
    print("O'k:"+str(ok*100/(ok+nok))+"%, Not O'k:"+str(nok*100/(ok+nok))+"%")