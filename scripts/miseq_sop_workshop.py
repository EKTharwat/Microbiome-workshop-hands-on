#!/usr/bin/env python3
"""Prepare or run the MiSeq SOP-style workshop on selected V3–V4 reads."""
import argparse
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile
from miseq import ROOT, install, resource


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['setup','run'])
    p.add_argument('--name',default='miseq-sop')
    p.add_argument('--processors',type=int,default=4)
    p.add_argument('--samples',nargs='+',help='Optional nonempty sample IDs for a small pilot')
    a=p.parse_args()
    if not a.name or Path(a.name).name!=a.name or a.name in ('.','..') or a.processors<1:p.error('Invalid run name or processor count')
    if a.command=='run' and a.samples is not None and len(set(a.samples))<3:p.error('The full SOP batch needs at least three samples for distance and PCoA steps')
    cmd=[sys.executable,str(ROOT/'scripts/workshop.py'),'setup','--name',a.name,'--processors',str(a.processors)]
    if a.samples:cmd+=['--samples',*a.samples]
    subprocess.run(cmd,check=True)
    work=ROOT/'runs'/a.name
    with (ROOT/'metadata/workshop_samples.tsv').open(newline='') as f:
        metadata={row['sample_name']:row for row in csv.DictReader(f,delimiter='\t')}
    selected=[line.split('\t',1)[0] for line in (work/'reads.files').read_text().splitlines()]
    with (work/'condition.design').open('w',newline='') as f:
        writer=csv.writer(f,delimiter='\t',lineterminator='\n');writer.writerow(['group','condition'])
        for sample in selected:writer.writerow([sample,metadata[sample]['condition']])
    archives=[('silva.bacteria.zip',['silva.bacteria/silva.bacteria.fasta']),('trainset9_032012.pds.zip',['trainset9_032012.pds.fasta','trainset9_032012.pds.tax'])]
    for archive,entries in archives:
        with zipfile.ZipFile(resource(archive)) as z:
            for entry in entries:
                with z.open(entry) as src,(work/Path(entry).name).open('wb') as dst:shutil.copyfileobj(src,dst)
    batch=(ROOT/'pipelines/miseq_sop_workshop.batch').read_text().replace('processors=4',f'processors={a.processors}',1)
    (work/'miseq_sop.batch').write_text(batch)
    print('SOP folder ready:',work,flush=True)
    if a.command=='setup':return
    exe_file=ROOT/'.tools/executable.txt'
    exe=Path(exe_file.read_text().strip()) if exe_file.exists() else install()
    if not exe.is_file():exe=install()
    env=dict(os.environ);env['PATH']=str(exe.parent)+os.pathsep+env.get('PATH','')
    with (work/'sop-console.log').open('w') as log:
        result=subprocess.run([str(exe),'miseq_sop.batch'],cwd=work,env=env,stdout=log,stderr=subprocess.STDOUT)
    text=(work/'sop-console.log').read_text(errors='replace')
    required=['final.opti_mcc.0.03.shared','final.asv.shared','final.taxonomy','final.tx.shared','final.opti_mcc.0.03.groups.ave-std.summary','final.opti_mcc.0.03.braycurtis.0.03.lt.ave.pcoa.axes']
    if result.returncode or '[ERROR]' in text or 'mothur > quit()' not in text or any(not (work/f).exists() for f in required):
        raise RuntimeError('SOP run incomplete; inspect '+str(work/'sop-console.log'))
    print('SOP run completed:',work)

if __name__=='__main__':
    try:main()
    except (OSError,ValueError,RuntimeError,StopIteration,zipfile.BadZipFile,subprocess.CalledProcessError) as exc:
        print('Error:',exc,file=sys.stderr);sys.exit(1)
