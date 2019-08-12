'''
Just to test that build is successful and built files are in place.
'''

from pathlib import Path
import sys

artifacts = (Path(__file__).parent / '../artifacts').absolute()
if sys.platform.startswith('linux'):
    assert (artifacts / 'lib/libReaktoro.so').is_file()
    assert (artifacts / 'include/Reaktoro/Reaktoro.hpp').is_file()
elif sys.platform == 'win32':
    assert (artifacts / 'lib/Reaktoro.lib').is_file()
    assert (artifacts / 'bin/Reaktoro.dll').is_file()
    assert (artifacts / 'include/Reaktoro/Reaktoro.hpp').is_file()
elif sys.platform == 'darwin':
    assert (artifacts / 'lib/libReaktoro.dylib').is_file()
    assert (artifacts / 'include/Reaktoro/Reaktoro.hpp').is_file()
else:
    assert False, f'Unsupported platform: {sys.platform}'

import reaktoro
assert hasattr(reaktoro, 'Reaction')
