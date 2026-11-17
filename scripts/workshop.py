#!/usr/bin/env python3
"""Prepare and run mothur on the bundled panel-filtered paired reads."""
import argparse,csv,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path
from miseq import ROOT,extract,install

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('command',choices=['setup','run'])
 p.add_argument('--name',default='workshop')
 p.add_argument('--processors',type=int,default=4)
 p.add_argument('--samples',nargs='+',help='Optional selected nonempty sample IDs')
 a=p.parse_args()
 if Path(a.name).name!=a.name or a.name in ('','.','..') or a.processors<1:p.error('invalid name or processors')
 archive=ROOT/'vendor/workshop_reads.zip'
 expected=json.loads((ROOT/'resources.json').read_text())['workshop_reads_sha256']
 if hashlib.sha256(archive.read_bytes()).hexdigest()!=expected:raise ValueError('Filtered data archive checksum mismatch')
 out=ROOT/'runs'/a.name
 out.mkdir(parents=True,exist_ok=False)
 extract(archive,out)
 with (out/'filtered_fastq/summary.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 selected=set(a.samples or [r['sample_id'] for r in rows if int(r['retained_pairs'])>0])
 valid={r['sample_id'] for r in rows if int(r['retained_pairs'])>0}
 if not selected or not selected<=valid:raise ValueError('Selected samples must exist and have retained pairs')
 with (out/'reads.files').open('w') as f:
  for row in rows:
   s=row['sample_id']
   if s in selected:f.write(f'{s}\tfiltered_fastq/{s}_R1.fastq.gz\tfiltered_fastq/{s}_R2.fastq.gz\n')
 shutil.copy2(ROOT/'references/panel/targets.16s.fasta',out/'targets.16s.fasta')
 with (ROOT/'metadata/panel/target_panel.tsv').open() as f:panel={r['target_id']:r['species'] for r in csv.DictReader(f,delimiter='\t')}
 with (out/'targets.taxonomy').open('w') as f:
  for line in (out/'targets.16s.fasta').read_text().splitlines():
   if line.startswith('>'):
    name=line[1:].split()[0];species=panel[name.split('__')[0]].replace(' ','_')
    f.write(f'{name}\tBacteria;{species};\n')
 (out/'pipeline.batch').write_text((ROOT/'pipeline.batch').read_text().replace('processors=4',f'processors={a.processors}',1))
 print(f'Prepared {out}',flush=True)
 if a.command=='setup':return
 exe_file=ROOT/'.tools/executable.txt'
 executable=Path(exe_file.read_text().strip()) if exe_file.exists() else install()
 if not executable.is_file():executable=install()
 env=dict(os.environ);env['PATH']=str(executable.parent)+os.pathsep+env.get('PATH','')
 with (out/'console.log').open('w') as log:result=subprocess.run([str(executable),'pipeline.batch'],cwd=out,stdout=log,stderr=subprocess.STDOUT,env=env)
 text=(out/'console.log').read_text(errors='replace')
 if result.returncode or '[ERROR]' in text or 'mothur > quit()' not in text:raise RuntimeError('mothur failed; inspect '+str(out/'console.log'))
 print('mothur complete:',out)
if __name__=='__main__':
 try:main()
 except (OSError,ValueError,RuntimeError,StopIteration) as e:print('Error:',e,file=sys.stderr);sys.exit(1)
