"""Verify opt-in ISA scope from real compile commands and compiler dependencies."""
import json
from pathlib import Path

build = Path('/tmp/flasheda-gpl-timing/build')
expected = {'nesterovBase.cpp', 'fft.cpp', 'fftsg.cpp', 'fftsg2d.cpp',
            'hpwl.cpp', 'wirelengthGradient.cpp', 'wirelengthSoa.cpp',
            'densityGradient.cpp'}
commands = json.loads((build / 'compile_commands.json').read_text())
seen = set()
for row in commands:
    if not row.get('output', '').startswith('src/gpl/CMakeFiles/gpl_lib.dir/'):
        continue
    name = Path(row['file']).name
    command = row['command']
    if name in expected:
        assert '-march=icelake-server' in command, name
        assert '-ffp-contract=off' in command and '-flto' in command, name
        deps = build / 'src/gpl/CMakeFiles/gpl_lib.dir/src' / (name + '.o.d')
        # Scope guard: kernels compiled with a different ISA must not instantiate
        # Eigen's ISA-dependent inline allocation/vectorization definitions.
        assert 'eigen3/' not in deps.read_text().lower(), name
        seen.add(name)
    else:
        assert '-march=icelake-server' not in command, name
assert seen == expected, (seen, expected)
print('Verified all eight kernel ISA flags, native ISA for other GPL sources, and no Eigen dependency in specialized kernels.')
