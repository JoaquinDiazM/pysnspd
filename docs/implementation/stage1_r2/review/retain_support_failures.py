"""Append observed support rejections; changing a grid must not erase a failure."""
from pathlib import Path
import argparse
import hashlib
import json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('assessments', type=Path, nargs='+')
    args = parser.parse_args()
    registry = Path(__file__).with_name('support_regression_points.json')
    history = json.loads(registry.read_text(encoding='utf-8'))
    known = {(point['amplitude'], point['gamma']) for point in history['points']}
    original_count = len(known)
    for path in args.assessments:
        assessment = json.loads(path.read_text(encoding='utf-8'))
        history['source_assessment_sha256'][path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        failures = assessment['dense_cell_support']['failures']
        failures += [row for row in assessment['exceptions'] if row.get('kind') == 'support']
        for row in failures:
            key = row['amplitude'], row['gamma']
            if key not in known:
                history['points'].append(dict(amplitude=key[0], gamma=key[1], first_failure_sources=[path.name]))
                known.add(key)
    registry.write_text(json.dumps(history, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(previous_points=original_count, retained_points=len(known),
                         registry_sha256=hashlib.sha256(registry.read_bytes()).hexdigest()), indent=2))


if __name__ == '__main__':
    main()
