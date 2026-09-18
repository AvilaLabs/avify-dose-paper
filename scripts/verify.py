#!/usr/bin/env python3
"""Reproduce historical numerical checks without OpenMC or external downloads."""
from pathlib import Path
import os, subprocess, sys, hashlib
ROOT=Path(__file__).resolve().parents[1]
for line in (ROOT/'MANIFEST.sha256').read_text().splitlines():
    digest, name = line.split('  ',1)
    if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
        raise SystemExit('Manifest mismatch: '+name)
env=os.environ.copy()
env['AUDIT_SCRATCH']=str(ROOT/'artifacts/audit')
env['W003ENV_PYTHON']=sys.executable
# Avoid silently using a local author archive; verify it only when explicitly supplied.
env.setdefault('W006_ARCHIVE',str(ROOT/'external/reference-archive-not-bundled.zip'))
subprocess.run([sys.executable,str(ROOT/'paper/manuscript/audit_numbers.py')],cwd=ROOT,env=env,check=True)
print('Public file hashes and recorded-result audit passed. Historical validation remains incomplete.')
