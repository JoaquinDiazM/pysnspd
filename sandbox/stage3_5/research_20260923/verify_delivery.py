"""Verify the documentary stage-3.5 delivery without writing or running physics.

The closed stage-3 delivery uses its archived command notebook. Its scientific
scope remains unchanged. Current research proposals are not physical admission.
"""
from pathlib import Path
import argparse
import json
import runpy


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT/'docs/implementation/stage3_5/research_20260923'
PREVIOUS = DATA/'previous_delivery_exact'
MANIFEST = DATA/'delivery_manifest.json'
OLD_MANIFEST = 'docs/implementation/stage3/closure_20260923/delivery_manifest.json'
OLD_VERIFIER = ROOT/'sandbox/stage3_spatial/closure_20260923/verify_delivery.py'
BASE = ROOT/'docs/implementation/stage3_5/parameter_inventory.json'
SELF = Path(__file__).resolve().relative_to(ROOT).as_posix()


def predecessors(old):
    """Preserve the earlier chain, then verify the 91-file stage-3 closure."""
    earlier_count = old['predecessors']()
    preserved = PREVIOUS/OLD_MANIFEST
    if old['sha'](preserved) != old['sha'](ROOT/OLD_MANIFEST):
        raise RuntimeError('The preserved stage-3 closure manifest differs from the frozen original')
    previous = old['read'](preserved)
    records = old['entries'](previous)
    if len(records) != 91:
        raise RuntimeError('Expected the original 91-file stage-3 closure delivery')
    overrides = {
        'docs/GEMINGA_COMMANDS.md': PREVIOUS/'docs/GEMINGA_COMMANDS.md',
        OLD_MANIFEST: preserved,
    }
    closure_count = old['verify'](records, overrides=overrides, require_size=True)
    for key in ('references', 'retained_references'):
        closure_count += old['verify'](previous.get(key, []), overrides=overrides)
    # Calling scope is intentional; old.main would compare the obsolete notebook.
    old['scope']()
    return earlier_count, closure_count


def require_zero(value, key, label):
    number = value.get(key)
    if isinstance(number, bool) or not isinstance(number, int) or number != 0:
        raise RuntimeError(f'{label}: {key} must be integer zero')


def require_acceptance(decisions, key):
    item = decisions.get(key)
    if not isinstance(item, dict) or item.get('accepted') is not True:
        raise RuntimeError('Missing explicit user acceptance: '+key)
    if not isinstance(item.get('answer'), str) or not item['answer'].strip():
        raise RuntimeError('Accepted decision lacks the recorded answer: '+key)
    return item


def sigma_range(value, allow_unset=False):
    """Read the declared nm scenario interval without assigning its bounds."""
    if isinstance(value, dict):
        if 'range' not in value and 'range_nm' not in value:
            raise RuntimeError('Gaussian scenario must explicitly declare or defer its range')
        if 'range' in value and 'range_nm' in value and value['range'] != value['range_nm']:
            raise RuntimeError('Ambiguous Gaussian scenario range')
        value = value.get('range', value.get('range_nm'))
    if value is None and allow_unset:
        return None
    if not isinstance(value, list) or len(value) != 2:
        raise RuntimeError('Gaussian scenario must declare its two range bounds in nm')
    return value


def transfer_scope(transfer, rows, proposal=None):
    """Keep the unapproved old family separate from source-based row proposals."""
    if not isinstance(transfer, dict) or transfer.get('accepted') is not False:
        raise RuntimeError('The historical transfer family was not approved')
    if transfer.get('proposed_values') is not None:
        raise RuntimeError('Superseded transfer-family values must remain historical')
    if not isinstance(transfer.get('status'), str) or not transfer['status'].strip():
        raise RuntimeError('Historical transfer-family research must retain its status')
    if proposal is None:
        # The documentary state before the new source review is still legitimate.
        for key in ('retention', 'Ein', 'photon_sigma', 'photon_position'):
            row = rows[key]
            status = row.get('status', '')
            pending = isinstance(status, str) and (
                'AWAITING' in status or 'PENDING' in status
                or status == 'USER_REQUESTED_FURTHER_PRIMARY_SOURCE_RESEARCH')
            if not pending or 'ACCEPTED' in status or 'ADOPTED' in status:
                raise RuntimeError('Source proposal missing for the transfer ledger: '+key)
            if row.get('value_or_proposal') is not None:
                raise RuntimeError('Unregistered transfer proposal: '+key)
        return False
    if not isinstance(proposal, dict):
        raise RuntimeError('Source-based transfer proposal must be an object')
    for key in ('retention', 'Ein_eV', 'gaussian_sigma_nm'):
        if key not in proposal or proposal[key] is None:
            raise RuntimeError('Source-based transfer proposal is incomplete: '+key)
    choice = proposal.get('user_choice')
    if not isinstance(choice, dict) or type(choice.get('accepted')) is not bool:
        raise RuntimeError('Conditional scenario requires an explicit acceptance boolean')
    accepted = choice['accepted']
    deferred = choice.get('status') == 'CHARACTERIZE_CASCADE_FIRST'
    if deferred and accepted:
        raise RuntimeError('Cascade-first decision cannot simultaneously accept a width scenario')
    if accepted and (not isinstance(choice.get('answer'), str) or not choice['answer'].strip()):
        raise RuntimeError('Accepted conditional scenario lacks the recorded user answer')
    if deferred and (not isinstance(choice.get('answer'), str) or not choice['answer'].strip()):
        raise RuntimeError('Cascade-first decision lacks the recorded user answer')
    sigma_status = ('USER_DEFERRED_PENDING_CASCADE' if deferred else
                    'USER_ACCEPTED_CONDITIONAL_SCENARIO' if accepted else
                    'CONDITIONAL_SCENARIO_PENDING_USER')
    expected = {
        'retention': 'SOURCE_REFERENCE_PROPOSED',
        'Ein': 'SOURCE_REFERENCE_PROPOSED',
        'photon_position': 'CONTROLLED_REFERENCE_PROPOSED',
        'photon_sigma': sigma_status,
    }
    for key, status in expected.items():
        if rows[key].get('status') != status:
            raise RuntimeError('Source proposal and parameter status disagree: '+key)
    proposed_range = sigma_range(proposal['gaussian_sigma_nm'], allow_unset=deferred)
    ledger_range = sigma_range(rows['photon_sigma'].get('value_or_proposal'), allow_unset=deferred)
    if deferred and (proposed_range is not None or ledger_range is not None):
        raise RuntimeError('Cascade-first decision must not activate a Gaussian width range')
    if ledger_range != proposed_range:
        raise RuntimeError('Gaussian scenario range differs between proposal and ledger')
    if 'physical_width_range_identified' in proposal and proposal['physical_width_range_identified'] is not False:
        raise RuntimeError('Source-based scenario does not identify a physical width range')
    return accepted


def research_scope(old):
    ledger = old['read'](DATA/'range_ledger.json')
    base = old['read'](BASE)
    decisions = old['read'](DATA/'user_decisions.json')
    if ledger.get('base_inventory_sha256') != old['sha'](BASE):
        raise RuntimeError('Range ledger is not tied to the unchanged base inventory')
    if ledger.get('user_decisions') != 'user_decisions.json':
        raise RuntimeError('Range ledger references a different decision record')
    items = ledger.get('parameters')
    base_items = base.get('parameters')
    if not isinstance(items, list) or not isinstance(base_items, list):
        raise RuntimeError('Parameter inventories require lists')
    ids = [item.get('id') for item in items]
    base_ids = [item.get('id') for item in base_items]
    if (ledger.get('count') != 127 or len(ids) != 127 or len(set(ids)) != 127
            or len(base_ids) != 127 or len(set(base_ids)) != 127
            or set(ids) != set(base_ids)):
        raise RuntimeError('Range ledger must cover exactly the 127 original parameter IDs')
    base_rows = {item['id']: item for item in base_items}
    rows = {item['id']: item for item in items}
    families = {item.get('family') for item in items}
    if ledger.get('families') != 15 or len(families) != 15:
        raise RuntimeError('Expected the 15 original parameter families')
    for key, row in rows.items():
        if row.get('family') != base_rows[key].get('family'):
            raise RuntimeError('Parameter family was changed: '+key)
        if 'physical_confidence_interval' not in row or row['physical_confidence_interval'] is not None:
            raise RuntimeError('Research ledger has silently admitted a physical interval: '+key)
    old['require_flags'](ledger, dict(independent_physical_confidence_box_admitted=False,
                                    code_parameters_changed=False), 'Research ledger')
    require_zero(ledger, 'new_physical_transients', 'Research ledger')
    require_acceptance(decisions, 'priority')
    time_origin = require_acceptance(decisions, 'time_origin')
    material = require_acceptance(decisions, 'material_reference')
    old['require_flags'](time_origin, dict(optical_cascade_delay_identified=False), 'Time origin')
    old['require_flags'](decisions, dict(no_automatic_physics_activation=True,
        long_computations_authorized_for_agent=False), 'Research decisions')
    if (rows['D']['value_or_proposal']['value'] != material.get('D_m2_s')
            or rows['sheet_resistance']['value_or_proposal']['value'] != material.get('sheet_resistance_ohm')):
        raise RuntimeError('Material ledger and accepted reference disagree')
    accepted = transfer_scope(decisions.get('controlled_transfer_family'), rows,
                              decisions.get('source_based_transfer_proposal'))
    return ledger, decisions, accepted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    old = runpy.run_path(str(OLD_VERIFIER))
    earlier_count, closure_count = predecessors(old)
    value = old['read'](MANIFEST)
    records = old['entries'](value)
    manifest_name = MANIFEST.relative_to(ROOT).as_posix()
    names = {item['path'] for item in records}
    if manifest_name in names:
        raise RuntimeError('The current delivery manifest must exclude itself')
    required = {SELF, 'docs/GEMINGA_COMMANDS.md',
                (DATA/'range_ledger.json').relative_to(ROOT).as_posix(),
                (DATA/'user_decisions.json').relative_to(ROOT).as_posix(),
                (PREVIOUS/'docs/GEMINGA_COMMANDS.md').relative_to(ROOT).as_posix(),
                (PREVIOUS/OLD_MANIFEST).relative_to(ROOT).as_posix()}
    if not required.issubset(names):
        raise RuntimeError('Delivery manifest omits contract files: '+', '.join(sorted(required-names)))
    count = old['verify'](records, require_size=True)
    reference_count = 0
    for key in ('references', 'retained_references'):
        reference_count += old['verify'](value.get(key, []))
    ledger, decisions, conditional_accepted = research_scope(old)
    for key in ('independent_physical_confidence_box_admitted', 'code_parameters_changed',
                'new_physics_activated', 'production_promotion'):
        if key in value:
            old['require_flags'](value, {key: False}, 'Delivery metadata')
    for key in ('new_physical_transients', 'new_simulations', 'physical_calculations_executed'):
        if key in value:
            require_zero(value, key, 'Delivery metadata')
    external = Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if external.exists() and old['sha'](external) != old['sha'](ROOT/'docs/GEMINGA_COMMANDS.md'):
        raise RuntimeError('External Geminga command notebook differs from the repository')
    print(json.dumps(dict(verified=True, files=count, references=reference_count,
        earlier_predecessor_records=earlier_count, stage3_closure_records=closure_count,
        manifest_sha256=old['sha'](MANIFEST), parameters=ledger['count'], families=ledger['families'],
        controlled_transfer_family_accepted=False,
        controlled_transfer_family_status=decisions['controlled_transfer_family'].get('status'),
        source_based_conditional_scenario_accepted=conditional_accepted,
        source_based_user_choice_status=decisions.get('source_based_transfer_proposal', {}).get('user_choice', {}).get('status'),
        physical_confidence_box_admitted=False, new_physical_transients=0,
        code_parameters_changed=False, production_promotion=False), indent=2))


if __name__ == '__main__':
    main()
