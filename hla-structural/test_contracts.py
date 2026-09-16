#!/usr/bin/env python3
"""Small substantive tests: exact markers, quality/mates and diploid optimization."""
import random, struct, subprocess, tempfile
from pathlib import Path
import numpy as np
from prepare_markers import kmers
from benchmark import fit_pair

ROOT=Path(__file__).resolve().parent

def main():
    rng=random.Random(17);s=''.join(rng.choices('ACGT',k=3000))
    keys=sorted(set(kmers(s)));assert keys
    rc=s.translate(str.maketrans('ACGT','TGCA'))[::-1]
    assert sorted(kmers(s))==sorted(kmers(rc))
    with tempfile.TemporaryDirectory() as tmp:
        p=Path(tmp);(p/'markers').write_text(''.join(str(x)+'\n' for x in keys))
        (p/'reads').write_text('@pair\n'+s+'\n+\n'+'I'*len(s)+'\n@pair\n'+rc+'\n+\n'+'I'*len(s)+'\n@bad\n'+s+'\n+\n'+'!'*len(s)+'\n')
        subprocess.run([str(ROOT/'count_markers'),'fastq',str(p/'markers'),str(p/'reads'),str(p/'counts')],check=True)
        assert set(struct.unpack('<'+'I'*len(keys),(p/'counts').read_bytes()))=={1}
        (p/'seq').write_text('>copies\n'+s+'N'+rc+'\n')
        subprocess.run([str(ROOT/'count_markers'),'fasta',str(p/'markers'),str(p/'seq'),str(p/'counts')],check=True)
        actual=np.fromfile(p/'counts',dtype='<u4')
        from collections import Counter
        expected=Counter(kmers(s+'N'+rc))
        assert actual.tolist()==[expected[k] for k in keys]
    X=np.array([[1,0,2],[0,1,1],[1,1,0]],dtype=float);y=X[0]+X[1]
    scores=fit_pair(X,y)
    assert np.unravel_index(np.argmin(scores),scores.shape)==(0,1)
    for i in range(3):
        for j in range(i,3):assert abs(scores[i,j]-np.sum((X[i]+X[j]-y)**2))<1e-8
    print('PASS: reverse complements, exact copy counts, Q20 filtering, mate deduplication, exhaustive pair scores')

if __name__=='__main__':main()
