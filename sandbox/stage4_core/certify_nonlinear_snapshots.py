"""Verify retained nonlinear snapshot bytes; no numerical solves."""
from pathlib import Path
import argparse,hashlib,json,time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--failed-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    start=time.monotonic();root=Path(__file__).resolve().parents[2];summary=read(a.root/'summary.json');identity=read(a.root/'identity.json')
    for name,digest in identity['sources'].items():
        if sha(root/name)!=digest:raise ValueError('Executed source changed: '+name)
    for name,digest in identity['inputs'].items():
        if sha(name)!=digest:raise ValueError('Frozen input changed: '+name)
    modes=[]
    for row in summary['records']:
        path=a.root/row['fields_path']
        if sha(path)!=row['fields_sha256']:raise ValueError('Output mode changed')
        modes.append(dict(n=row['n'],path=row['fields_path'],sha256=row['fields_sha256'],bytes=path.stat().st_size))
    failed=[]
    for path in sorted((a.failed_root/'modes').glob('*.json')):
        row=read(path);fields=a.failed_root/row['fields_path']
        if sha(fields)!=row['fields_sha256']:raise ValueError('Preserved initial checkpoint changed')
        failed.append(dict(n=row['n'],path=row['fields_path'],sha256=row['fields_sha256'],bytes=fields.stat().st_size))
    data=dict(status='NONLINEAR_SNAPSHOTS_CERTIFIED',runtime_seconds=time.monotonic()-start,
        full_root=str(a.root),frequency_modes=modes,total_mode_bytes=sum(x['bytes'] for x in modes),
        compact_files={name:sha(a.root/name) for name in ('summary.json','identity.json','executed_plan.json','nonlinear_comparison_fields.npz')},
        original_failed_root=str(a.failed_root),original_failure=read(a.failed_root/'failure.json'),
        original_identity_sha256=sha(a.failed_root/'identity.json'),original_preserved_modes=failed,
        source_hashes_verified=len(identity['sources']),input_hashes_verified=len(identity['inputs']),
        certificate_source_sha256=sha(Path(__file__)),spectral_tolerance_changed=False,physical_equations_changed=False)
    with a.output.open('x') as stream:json.dump(data,stream,indent=2);stream.write('\n')
    print(json.dumps({key:value for key,value in data.items() if key not in ('frequency_modes','original_preserved_modes')},indent=2))


if __name__=='__main__':main()
