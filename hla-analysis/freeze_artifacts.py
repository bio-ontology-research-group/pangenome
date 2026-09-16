#!/usr/bin/env python3
"""Write final local artifact and software hashes after the workflow/audit finish."""
from pathlib import Path
import hashlib
import json
import platform
import sys
import Bio
import numpy
import scipy
import matplotlib
import openpyxl


def main():
    base=Path(__file__).resolve().parent;files=[]
    for p in sorted(base.rglob('*')):
        if not p.is_file() or '__pycache__' in p.parts or p.name=='artifact_manifest.json':continue
        files.append(dict(path=p.relative_to(base).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    data=dict(files=files,python=sys.version,platform=platform.platform(),biopython=Bio.__version__,numpy=numpy.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,openpyxl=openpyxl.__version__)
    (base/'results/artifact_manifest.json').write_text(json.dumps(data,indent=2)+'\n')
    print(len(files),'local artifacts hashed')


if __name__=='__main__':main()
