"""Read-only verification of the scoped 3.5 closure and non-photon stage-4 preparation.

Archived front doors are used for the previous delivery. This verifier never
generates manifests, edits decisions, launches physics or promotes production.
"""
from pathlib import Path
import argparse
import json
import runpy
import sys


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT/'docs/implementation/stage3_5/assessment_r2_20260923'
PREVIOUS = DATA/'previous_delivery_exact'
MANIFEST = DATA/'delivery_manifest.json'
R1_MANIFEST = 'docs/implementation/stage3_5/research_20260923/delivery_manifest.json'
R1_LEDGER = ROOT/'docs/implementation/stage3_5/research_20260923/range_ledger.json'
R1_VERIFIER = ROOT/'sandbox/stage3_5/research_20260923/verify_delivery.py'
STAGE3_VERIFIER = ROOT/'sandbox/stage3_spatial/closure_20260923/verify_delivery.py'
BASE = ROOT/'docs/implementation/stage3_5/parameter_inventory.json'
SELF = Path(__file__).resolve().relative_to(ROOT).as_posix()
USER_ANSWER = 'Cerrar investigación 3.5 y preparar etapa 4 sin fotón'
EXPECTED_CLOSURE = dict(research_closed=True, physical_domain_complete=False,
    photon_transfer_admitted=False, stage4_prepared=True, stage4_started=False,
    physics_activated=False, production_promotion=False)


def predecessors(stage3, r1):
    """Use the existing chain exactly once, then check the 36-file r1 payload."""
    earlier_count, stage3_count = r1['predecessors'](stage3)
    r1['research_scope'](stage3)
    preserved = PREVIOUS/R1_MANIFEST
    if stage3['sha'](preserved) != stage3['sha'](ROOT/R1_MANIFEST):
        raise RuntimeError('Preserved r1 manifest differs from the frozen original')
    previous = stage3['read'](preserved)
    records = stage3['entries'](previous)
    if len(records) != 36:
        raise RuntimeError('Expected the original 36-file r1 delivery')
    overrides = {name: PREVIOUS/name for name in (
        'docs/GEMINGA_COMMANDS.md', 'docs/implementation/stage3_5/CURRENT.md', R1_MANIFEST)}
    count = stage3['verify'](records, overrides=overrides, require_size=True)
    for key in ('references', 'retained_references'):
        count += stage3['verify'](previous.get(key, []), overrides=overrides, require_size=True)
    return earlier_count, stage3_count, count


def require_zero(value, key, label):
    number = value.get(key)
    if type(number) is not int or number != 0:
        raise RuntimeError(f'{label}: {key} must be integer zero')


def validate_scope(closure, domain, base, require_flags):
    require_flags(closure, EXPECTED_CLOSURE, 'Stage3.5 scoped closure')
    require_zero(closure, 'new_physical_transients', 'Stage3.5 scoped closure')
    if closure.get('user_answer') != USER_ANSWER:
        raise RuntimeError('The closure must retain the explicit recorded user answer')
    if closure.get('stage4_scope') != 'NON_PHOTON_DIAGNOSTICS':
        raise RuntimeError('Stage4 preparation must be restricted to non-photon diagnostics')
    require_flags(domain, dict(physical_confidence_box_admitted=False,
        code_parameters_changed=False, all_solver_parameters_unchanged=True), 'Use domain')
    require_zero(domain, 'new_transients', 'Use domain')
    if 'photon_width_nm' not in domain or domain['photon_width_nm'] is not None:
        raise RuntimeError('The physical photon width must remain explicitly open')
    items = domain.get('parameters')
    base_items = base.get('parameters')
    if not isinstance(items, list) or not isinstance(base_items, list):
        raise RuntimeError('Use domain and base inventory require parameter lists')
    ids = [item.get('id') for item in items]
    base_ids = [item.get('id') for item in base_items]
    if (domain.get('count') != 127 or len(ids) != 127 or len(set(ids)) != 127
            or len(base_ids) != 127 or len(set(base_ids)) != 127
            or set(ids) != set(base_ids)):
        raise RuntimeError('Use domain must cover exactly the 127 original parameter IDs')
    families = {item.get('family') for item in items}
    if domain.get('families') != 15 or len(families) != 15:
        raise RuntimeError('Use domain must retain the 15 original parameter families')
    originals = {item['id']: item for item in base_items}
    rows = {item['id']: item for item in items}
    for key, row in rows.items():
        if row.get('family') != originals[key].get('family'):
            raise RuntimeError('Parameter family changed: '+key)
        if 'physical_confidence_interval' not in row or row['physical_confidence_interval'] is not None:
            raise RuntimeError('Use domain silently admits a physical interval: '+key)
        if 'implementation_activation' in row and row['implementation_activation'] is not False:
            raise RuntimeError('Research activated a solver parameter: '+key)
    width = rows['photon_sigma'].get('value_or_proposal')
    if not isinstance(width, dict) or 'value' not in width or 'range' not in width:
        raise RuntimeError('Photon-width row must explicitly retain value:null and range:null')
    if width['value'] is not None or width['range'] is not None:
        raise RuntimeError('Physical photon width was selected despite the cascade-first decision')
    return rows


def scope(stage3):
    closure = stage3['read'](DATA/'closure_decision.json')
    domain = stage3['read'](DATA/'use_domain.json')
    base = stage3['read'](BASE)
    if domain.get('base_inventory_sha256') != stage3['sha'](BASE):
        raise RuntimeError('Use domain is not tied to the unchanged base inventory')
    if domain.get('r1_ledger_sha256') != stage3['sha'](R1_LEDGER):
        raise RuntimeError('Use domain is not tied to the preserved r1 ledger')
    validate_scope(closure, domain, base, stage3['require_flags'])
    return closure, domain


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    stage3 = runpy.run_path(str(STAGE3_VERIFIER))
    r1 = runpy.run_path(str(R1_VERIFIER))
    earlier_count, stage3_count, r1_count = predecessors(stage3, r1)
    manifest = stage3['read'](MANIFEST)
    records = stage3['entries'](manifest)
    names = {item['path'] for item in records}
    if MANIFEST.relative_to(ROOT).as_posix() in names:
        raise RuntimeError('The current manifest must exclude itself')
    required = {SELF, 'docs/GEMINGA_COMMANDS.md',
        'docs/implementation/stage3_5/CURRENT.md',
        'docs/implementation/stage4/entry_contract.json',
        'docs/implementation/stage4/README.md',
        (DATA/'closure_decision.json').relative_to(ROOT).as_posix(),
        (DATA/'use_domain.json').relative_to(ROOT).as_posix()}
    required.update((PREVIOUS/name).relative_to(ROOT).as_posix() for name in (
        'docs/GEMINGA_COMMANDS.md', 'docs/implementation/stage3_5/CURRENT.md', R1_MANIFEST))
    if not required.issubset(names):
        raise RuntimeError('Delivery omits contract files: '+', '.join(sorted(required-names)))
    count = stage3['verify'](records, require_size=True)
    reference_count = 0
    for key in ('references', 'retained_references'):
        reference_count += stage3['verify'](manifest.get(key, []), require_size=True)
    closure, domain = scope(stage3)
    for key, flag in EXPECTED_CLOSURE.items():
        if key in manifest:
            stage3['require_flags'](manifest, {key: flag}, 'Delivery metadata')
    if 'new_physical_transients' in manifest:
        require_zero(manifest, 'new_physical_transients', 'Delivery metadata')
    external = Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if sys.platform.startswith('linux') and external.exists():
        if stage3['sha'](external) != stage3['sha'](ROOT/'docs/GEMINGA_COMMANDS.md'):
            raise RuntimeError('External Geminga command notebook differs from the repository')
    print(json.dumps(dict(verified=True, files=count, references=reference_count,
        predecessor_records=dict(earlier=earlier_count, stage3_closure=stage3_count, research_r1=r1_count),
        manifest_sha256=stage3['sha'](MANIFEST), parameters=domain['count'], families=domain['families'],
        research_closed=closure['research_closed'], stage4_prepared=True, stage4_started=False,
        stage4_scope=closure['stage4_scope'], photon_width_nm=None,
        physical_domain_complete=False, photon_transfer_admitted=False,
        new_physical_transients=0, physics_activated=False, production_promotion=False), indent=2))


if __name__ == '__main__':
    main()
