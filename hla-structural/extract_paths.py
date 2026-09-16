#!/usr/bin/env python3
"""Extract actual graph traversals and verify their sequence against input assemblies."""
import csv, gzip, hashlib, json, re
from pathlib import Path
from collections import defaultdict

BASE=Path('/home/leechuck/hla/mhc_all')
GRAPH=Path('/home/asianhla/data/upload/HLA/mhc_graph/MHC.full.gfa.gz')
COMP=str.maketrans('ACGTN','TGCAN')

def fasta(path):
    out={}; name=None
    for line in open(path):
        if line.startswith('>'): name=line[1:].split()[0];out[name]=[]
        else: out[name].append(line.strip().upper())
    return {k:''.join(v) for k,v in out.items()}

def main():
    rows=list(csv.DictReader(open('intervals.tsv'),delimiter='\t'))
    req=defaultdict(list)
    for r in rows: req[r['graph_path']].append(r)
    nodes={}; walks={}; types=defaultdict(int)
    with gzip.open(GRAPH,'rt') as f:
        for line in f:
            x=line.rstrip().split('\t');types[x[0]]+=1
            if x[0]=='S': nodes[x[1]]=x[2].upper()
            elif x[0]=='W':
                name='#'.join(x[1:4])
                if name in req:
                    walks[name]=(int(x[4]) if x[4]!='*' else 0,re.findall(r'([><])([^><]+)',x[6]))
            elif x[0]=='P' and x[1] in req:
                walks[x[1]]=(0,[('>' if s[-1]=='+' else '<',s[:-1]) for s in x[2].split(',')])
    print('graph',dict(types),'requested paths',len(req),'found',len(walks),flush=True)
    result=[]; cache={}; graph_records=[]
    with gzip.open('locus_sequences.fa.gz','wt') as fa:
        for path, requests in req.items():
            if path not in walks: raise ValueError('Missing graph path '+path)
            offset, walk=walks[path]
            positions=[];pos=offset
            for orient,node in walk:
                seq=nodes[node];positions.append((pos,pos+len(seq),orient,node));pos+=len(seq)
            for r in requests:
                lo,hi=int(r['start0']),int(r['end0']);parts=[];steps=[]
                for a,b,orient,node in positions:
                    if a>=hi:break
                    if b<=lo:continue
                    seq=nodes[node]
                    if orient=='<':seq=seq.translate(COMP)[::-1]
                    start=max(lo,a)-a;end=min(hi,b)-a
                    parts.append(seq[start:end]);steps.append(dict(node=node,orientation=orient,start_offset=start,end_offset=end,length=end-start))
                seq=''.join(parts)
                entry=r['hap_id'].replace('#','_')
                if entry not in cache: cache={entry:fasta(BASE/(entry+'.mhc.fa'))}
                source=cache[entry][r['contig']][lo:hi]
                r['sequence_matches_assembly']=int(seq==source and len(seq)==hi-lo)
                r['N_bases']=seq.count('N');r['other_ambiguous_bases']=sum(c not in 'ACGTN' for c in seq)
                r['sequence_sha256']=hashlib.sha256(seq.encode()).hexdigest()
                r['sequence_id']=r['locus']+'__'+r['hap_id']
                r['eligible']=int(r['sequence_matches_assembly'] and not r['N_bases'] and not r['other_ambiguous_bases'])
                r['reads_available']=int(Path('/home/asianhla/data/upload/1000G_MHC/cram',r['donor']+'.mhc.cram').is_file())
                fa.write('>'+r['sequence_id']+'\n'+seq+'\n')
                graph_records.append(dict(sequence_id=r['sequence_id'],steps=steps))
                result.append(r)
    with open('catalogue.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
    with gzip.open('graph_paths.jsonl.gz','wt') as f:
        for r in graph_records:f.write(json.dumps(r)+'\n')
    Path('extraction_summary.json').write_text(json.dumps(dict(graph=str(GRAPH),record_counts=dict(types),intervals=len(result),eligible=sum(r['eligible'] for r in result),sequence_mismatches=sum(not r['sequence_matches_assembly'] for r in result)),indent=2)+'\n')
    print(Path('extraction_summary.json').read_text(),flush=True)

if __name__=='__main__':main()
