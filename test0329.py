#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import re
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.normalizer

def main(argv):
    trans = argv[1]
    anno = argv[2]
    output = argv[3]
    hp = hgvs.parser.Parser()
    hdp = hgvs.dataproviders.uta.connect()
    am = hgvs.assemblymapper.AssemblyMapper(hdp,assembly_name='GRCh37', alt_aln_method='splign',replace_reference=True)
    hn = hgvs.normalizer.Normalizer(hdp)
    my_d = {}
    my_strand={}
    with open(trans,'r') as NM :
        for i,hang in enumerate(NM,1):
            hang = hang.rstrip()
            t_version = hang.split('\t')
            my_d[t_version[0]] = t_version[1]
            my_strand[t_version[0]] = t_version[2]

    with open(anno,'r') as IN,\
        open(output,'w') as OUT :
        for n,line in enumerate(IN,1):
            line = line.rstrip()
            if line.startswith("Chr"):
                OUT.write("Chr_c\tStart_c\tEnd_c\tStrand\tTranscript_c\tHGVS_C_c\tHGVS_P_c\t"+str(line) +"\n")
            else:
                content = line.split('\t')
                transcrip=content[6]
                c_change=content[7]
                if transcrip in my_d :
                    Strand=my_strand[transcrip]
                    transcrip = my_d[transcrip]
                    #print transcrip
                    try:
                        c1 = hp.parse_hgvs_variant(transcrip + ':' + c_change) 
                        c1n = hn.normalize(c1)
                        g = am.c_to_g(c1)
                        p = am.c_to_p(c1)
                        pep = re.sub('[()]','',str(p))
                        transcrip,c_change = str(c1n).split(':')
                        np,p_change = str(pep).split(':')
                        start_pos=str(g.posedit.pos.start)
                        end_pos=str(g.posedit.pos.end)
                    except Exception,e:
                        print e.message
                        c_change="NA"
                        p_change="NA"
                        start_pos="NA"
                        end_pos=="NA"
                        continue
                    old='\t'.join(content)
                    OUT.write(content[0]+"\t"+str(start_pos)+"\t"+str(end_pos)+"\t"+Strand+"\t"+transcrip+"\t"+c_change+"\t"+p_change+"\t"+old+"\n")
                            
if __name__ == "__main__":
     main(sys.argv)
