#!/usr/bin/env python3
"""Recompute analysis from frozen source data and archived remote measurements.

Does not silently re-download references, rebuild the graph, or resubmit read jobs.
Remote recomputation instructions and scripts are in REPRODUCIBILITY.md.
"""
from pathlib import Path
import os
import subprocess
import sys
import time
import json

BASE=Path(__file__).resolve().parent
STEPS=['sequence_catalogue.py','coding_candidates.py','benchmark.py',
       'within_type_catalogue.py','experimental_discordances.py',
       'benchmark_uncertainty.py','bundle_increment.py','freeze_panel.py',
       'assess_contrasts.py','prepare_read_markers.py','prepare_flank_markers.py',
       'summarize_read_support.py','independent_assembly_checks.py',
       'published_sequence_check.py','rccx_catalogue.py','rccx_report.py',
       'summarize_additional_validation.py','check_results.py','final_report.py']


def main():
    env=dict(os.environ,OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
    logs=BASE/'results/reproduction';logs.mkdir(exist_ok=True)
    runs=[]
    for script in STEPS:
        start=time.monotonic();print('Running',script,flush=True)
        with (logs/(script+'.log')).open('w') as f:
            result=subprocess.run([sys.executable,str(BASE/script)],env=env,stdout=f,stderr=subprocess.STDOUT)
        runs.append(dict(script=script,exit_code=result.returncode,elapsed_seconds=round(time.monotonic()-start,3)))
        (logs/'run.json').write_text(json.dumps(runs,indent=2)+'\n')
        if result.returncode:raise SystemExit(f'Failed {script}; inspect {logs}')
    print('All analysis stages completed successfully',flush=True)


if __name__=='__main__':main()
