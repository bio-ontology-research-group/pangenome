#!/usr/bin/env python3
"""Standalone scientific figure; all attempted donors stay in denominators."""
import csv,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent

def main():
    rows=list(csv.DictReader(open(ROOT/'results/split_summary.tsv'),delimiter='\t'))
    calibrated=json.loads((ROOT/'results/calibrated_summary.json').read_text())
    rows.append(dict(locus='RCCX',method='calibrated',evaluation_set='additional_validation',**calibrated))
    fig,axes=plt.subplots(1,2,figsize=(12,5.6),sharey=True)
    methods=['path_multiplicity','binary_paths','calibrated','assembly_profile_oracle']
    labels=['Uncalibrated path dosage','Binary path markers','Pilot-calibrated CN + paths','Perfect assembly profile']
    colors=['#2471a3','#b9770e','#27845b','#68737d']
    for ax,locus in zip(axes,['RCCX','DRB']):
        for j,(m,label,color) in enumerate(zip(methods,labels,colors)):
            if locus=='DRB' and m=='calibrated':continue
            r=next(r for r in rows if r['locus']==locus and r['method']==m and r['evaluation_set']=='additional_validation')
            values=[int(r[k])/int(r['n']) for k in ['cn_correct','pair_correct']]
            ax.bar(np.arange(2)+(j-1.5)*.2,values,width=.18,color=color,label=label)
            for x,y,k in zip(np.arange(2)+(j-1.5)*.2,values,['cn_correct','pair_correct']):
                ax.text(x,y+.015,f"{r[k]}/{r['n']}",ha='center',va='bottom',fontsize=8,rotation=90)
        ax.set_xticks([0,1],['Total copy number','Structural pair']);ax.set_title(locus)
        ax.spines[['top','right']].set_visible(False);ax.set_ylim(0,1.18)
    axes[0].set_ylabel('Agreement with held-out assembly labels')
    fig.suptitle('Structural typing: 88 additional donors, family-excluded references',fontsize=13)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False)
    fig.text(.5,.13,'Calibrated reference depth also recovers RCCX total copies in 88/88; it does not resolve arrangement.',ha='center',fontsize=9)
    fig.text(.5,.10,'Perfect assembly profile is a diagnostic. Dosage calibration applies only to RCCX.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.17,1,.94]);fig.savefig(ROOT/'results/benchmark.png',dpi=180);fig.savefig(ROOT/'results/benchmark.pdf');plt.close(fig)

if __name__=='__main__':main()
