#!/usr/bin/env python3
"""Check antigen-recognition-exon compatibility, without rewriting Sanger truth."""
from collections import defaultdict, Counter
import json
import re
from sequence_catalogue import BASE, OUT, write_table
from check_results import table


def main():
    allele_groups={}; prefix_groups=defaultdict(set)
    for line in (BASE/'source/imgt/wmda/hla_nom_g.txt').read_text().splitlines():
        if line.startswith('#') or not line: continue
        locus,alleles,group=line.split(';')
        members=alleles.split('/')
        gid=locus+(group or members[0])
        for a in members:
            allele_groups[locus+a]=gid
            parts=re.sub('[A-Z]+$','',a).split(':')
            for n in range(1,len(parts)+1):
                prefix_groups[locus+':'.join(parts[:n])].add(gid)
    seq={r['name']:r for r in table(OUT/'sequence_catalogue.tsv')}
    out=[]
    for r in table(OUT/'experimental_genotypes.tsv'):
        if r['method']!='current_exact_cds' or r['match']!='0': continue
        locus=r['gene'].removeprefix('HLA-')+'*'
        observed=[]; annotated=[]
        for truth in r['truth'].split(';'):
            observed.append(set().union(*(prefix_groups.get(locus+re.sub('[A-Z]+$','',a),set()) for a in truth.split('/'))))
        source=[]
        for hap in ['1','2']:
            entry=seq[r['donor']+'#'+hap+'#'+r['gene']]
            source.append(entry['old_consensus'])
            annotated.append({allele_groups[a] for a in entry['exact_cds_alleles'].split(';') if a in allele_groups})
        possible=bool((observed[0]&annotated[0] and observed[1]&annotated[1]) or (observed[0]&annotated[1] and observed[1]&annotated[0]))
        state='compatible_with_antigen_exon_resolution' if possible else 'unresolved_antigen_exon_discordance'
        out.append(dict(donor=r['donor'],gene=r['gene'],experimental=r['truth'],current_exact_cds=r['prediction'],old_annotation=';'.join(source),assessment=state,experimental_possible_g_groups=json.dumps([sorted(s) for s in observed]),assembly_g_groups=json.dumps([sorted(s) for s in annotated]),limitation='G-group compatibility is possible exon identity under expanded historical ambiguity, not proof of assay error or whole-gene agreement; unresolved cases require read/assembly inspection'))
    write_table(OUT/'experimental_discordances.tsv',out)
    summary=dict(Counter(r['assessment'] for r in out))
    (OUT/'experimental_discordances_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary)


if __name__=='__main__':
    main()
