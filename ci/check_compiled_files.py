'''
Just to test that build is successful and built files are in place.
'''

from pathlib import Path
import os
import sys
from distutils.sysconfig import get_config_var


ext_suffix = get_config_var('EXT_SUFFIX')
pdb_suffix = '.'.join(ext_suffix.split('.')[:-1] + ['pdb'])

artifacts = (Path(__file__).parent / '../artifacts').absolute()

if sys.platform.startswith('linux'):
    if (artifacts / 'lib64').is_dir():
        assert (artifacts / 'lib64/libReaktoro.so').is_file()
    elif (artifacts / 'lib').is_dir():
        assert (artifacts / 'lib/libReaktoro.so').is_file()
    else:
        raise RuntimeError(f"Both {artifacts / 'lib64/libReaktoro.so'} " 
            f"and {artifacts / 'lib/libReaktoro.so'} are not found.")
    assert (artifacts / 'include/Reaktoro/Reaktoro.hpp').is_file()

elif sys.platform == 'win32':
    assert (artifacts / 'lib/Reaktoro.lib').is_file()
    assert (artifacts / 'bin/Reaktoro.dll').is_file()
    assert (artifacts / 'include/Reaktoro/Reaktoro.hpp').is_file()
    assert (artifacts / f'python/Lib/site-packages/reaktoro/PyReaktoro{ext_suffix}').is_file()
    if os.environ.get('CONFIG', '').lower() == 'debug':
        assert (artifacts / 'bin/Reaktoro.pdb').is_file()
        assert (artifacts / f'python/Lib/site-packages/reaktoro/PyReaktoro{pdb_suffix}').is_file()

elif sys.platform == 'darwin':
    assert (artifacts / 'lib/libReaktoro.dylib').is_file()
    assert (artifacts / 'include/Reaktoro/Reaktoro.hpp').is_file()

else:
    assert False, f'Unsupported platform: {sys.platform}'

import reaktoro
assert hasattr(reaktoro, 'Reaction')
