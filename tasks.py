from invoke import task
from pathlib import Path
import os
import os.path
import shutil
import sys


def strip_and_join(s):
    return ' '.join(line.strip() for line in s.splitlines() if line.strip() != '')


def echo(c, msg):
    from colorama.ansi import Fore, Style
    if c.config.run.echo:
        print(f"{Fore.WHITE}{Style.BRIGHT}{msg}{Style.RESET_ALL}")


def remove_directory(c, path):
    if path.is_dir():
        echo(c, f"Removing {path}")
        shutil.rmtree(path)
    else:
        echo(c, f"Not removing {path} (not a directory)")


if sys.platform.startswith('win'):

    @task
    def msvc(c, clean=False):
        """
        Generates a Visual Studio project at the "build" directory, ready to
        be used for development.
        Assumes that the environment is already configured using:
            conda devenv
            activate reaktoro
        """
        root_dir = Path(__file__).parent
        build_dir = root_dir / "build/msvc"
        artifacts_dir = root_dir / "artifacts"
        if clean:
            remove_directory(c, build_dir)
            remove_directory(c, artifacts_dir)
        build_dir.mkdir(parents=True, exist_ok=True)
        os.chdir(build_dir)

        relative_root_dir = Path(os.path.relpath(root_dir, Path.cwd()))
        relative_artifacts_dir = Path(os.path.relpath(artifacts_dir, Path.cwd()))

        CMAKE_GENERATOR = "Visual Studio 14 2015"
        CMAKE_ARCH = "x64"
        CONFIG = "Release"
        CMAKE_INCLUDE_PATH = f"{os.environ['CONDA_PREFIX']}\\Library\\include"

        # `PYTHON_INSTALL_PREFIX` is configured to `artifacts_dir / 'python'`
        # so that it won't "pollute" the Python environment when in develop
        # mode.

        c.run(strip_and_join(f"""
            cmake -G "{CMAKE_GENERATOR}"
            -A "{CMAKE_ARCH}"
            -DBUILD_ALL=ON
            -DPYTHON_INSTALL_PREFIX="{(relative_artifacts_dir / 'python').as_posix()}"
            -DCMAKE_BUILD_TYPE={CONFIG}
            -DCMAKE_INCLUDE_PATH="{CMAKE_INCLUDE_PATH}"
            -DBOOST_INCLUDE_DIR="{CMAKE_INCLUDE_PATH}"
            -DCMAKE_INSTALL_PREFIX="{relative_artifacts_dir.as_posix()}"
            -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
            "-DTHIRDPARTY_COMMON_ARGS=-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON"
            "{str(relative_root_dir)}"
        """))
