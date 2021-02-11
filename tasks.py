from invoke import task
from invoke.exceptions import Exit
from pathlib import Path
from typing import Optional
import os
import shutil
import sys


def _get_vcvars_paths():
    template = r"%PROGRAMFILES(X86)%\Microsoft Visual Studio\2019\{edition}\VC\Auxiliary\Build\vcvarsall.bat"
    template = os.path.expandvars(template)
    editions = ('BuildTools', 'Professional', 'WDExpress', 'Community', 'Enterprise')
    return tuple(Path(template.format(edition=edition)) for edition in editions)


def strip_and_join(s: str):
    return ' '.join(line.strip() for line in s.splitlines() if line.strip() != '')


def echo(c, msg: str):
    from colorama.ansi import Fore, Style
    if c.config.run.echo:
        print(f"{Fore.WHITE}{Style.BRIGHT}{msg}{Style.RESET_ALL}")


def remove_directory(c, path: Path):
    if path.is_dir():
        echo(c, f"Removing {path}")
        shutil.rmtree(path)
    else:
        echo(c, f"Not removing {path} (not a directory)")


def _get_and_prepare_build_and_artifacts_dir(
    c,
    clean: bool=False,
    build_subdirectory: str="build",
    ) -> Path:
    '''
    Returns build directory where `cmake` shall be called from. Creates it and
    possibly removes its contents (and artifacts_dir contents) if `clean=True`
    is passed.
    '''
    root_dir = Path(__file__).parent
    build_dir = root_dir / build_subdirectory
    artifacts_dir = root_dir / "artifacts"
    if clean:
        remove_directory(c, build_dir)
        remove_directory(c, artifacts_dir)
    build_dir.mkdir(parents=True, exist_ok=not clean)
    artifacts_dir.mkdir(parents=True, exist_ok=not clean)
    return build_dir, artifacts_dir


def _get_cmake_command(
    c,
    build_dir: Path,
    artifacts_dir: Path,
    cmake_generator: str,
    cmake_arch: Optional[str]=None,
    config: str='Release',
    verbose=False,
    ):
    '''
    :param build_dir: Directory from where cmake will be called.
    :param artifacts_dir: Directory where binaries and python modules will be installed.
    '''
    root_dir = Path(__file__).parent
    relative_root_dir = Path(os.path.relpath(root_dir, build_dir))
    relative_artifacts_dir = Path(os.path.relpath(artifacts_dir, build_dir))

    if sys.platform.startswith('win'):
        cmake_include_path = f"{os.environ['CONDA_PREFIX']}\\Library\\include"
    else:
        cmake_include_path = f"{os.environ['CONDA_PREFIX']}\\include"

    # `PYTHON_INSTALL_PREFIX` is configured to `artifacts_dir / 'python'`
    # so that it won't "pollute" the Python environment when in develop
    # mode.

    return strip_and_join(f"""
        cmake
            {f'-G "{cmake_generator}"' if cmake_generator is not None else ""}
            {f'-A "{cmake_arch}"' if cmake_arch is not None else ""}
            -DPYBIND11_PYTHON_VERSION={os.environ.get("PY_VER", "3.7")}
            -DREAKTORO_BUILD_ALL=ON            
            -DREAKTORO_PYTHON_INSTALL_PREFIX="{(artifacts_dir / 'python').as_posix()}"
            -DCMAKE_BUILD_TYPE={config}
            -DCMAKE_INCLUDE_PATH="{cmake_include_path}"
            -DCMAKE_INSTALL_PREFIX="{relative_artifacts_dir.as_posix()}"
            {f'-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON' if verbose else ''}
            {f'-DREAKTORO_USE_OPENLIBM:BOOL=ON' if os.environ.get('REAKTORO_USE_OPENLIBM', '0').strip() != '0' else ''}
            {f'"-DREAKTORO_THIRDPARTY_EXTRA_BUILD_ARGS=-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON"' if verbose else ''}
            "{str(relative_root_dir)}"
    """)


if sys.platform.startswith('win'):
    @task
    def msvc(c, clean=False, config='Release'):
        """
        Generates a Visual Studio project at the "build/msvc" directory.
        Assumes that the environment is already configured using:
            conda devenv
            activate reaktoro
        """
        build_dir, artifacts_dir = _get_and_prepare_build_and_artifacts_dir(
            c,
            clean=clean,
            build_subdirectory="build/msvc",
        )
        cmake_command = _get_cmake_command(
            c,
            build_dir=build_dir,
            artifacts_dir=artifacts_dir,
            cmake_generator="Visual Studio 16 2019",
            cmake_arch="x64",
            config=config,
        )
        os.chdir(build_dir)
        c.run(cmake_command)


@task
def compile(c, clean=False, config='Release', number_of_jobs=-1, verbose=False):
    """
    Compiles Reaktoro by running CMake and building with `ninja`.
    Assumes that the environment is already configured using:
        conda devenv
        [source] activate reaktoro
    """
    build_dir, artifacts_dir = _get_and_prepare_build_and_artifacts_dir(
        c,
        clean=clean,
        build_subdirectory="build",
    )

    cmake_command = _get_cmake_command(
        c,
        build_dir=build_dir,
        artifacts_dir=artifacts_dir,
        cmake_generator="Ninja",
        config=config,
        verbose=verbose,
    )
    build_command = strip_and_join(f"""
        cmake
            --build .
            --target install
            --config {config}
            --
                {f"-j {number_of_jobs}" if number_of_jobs > 0 else ""}
    """)

    commands = [cmake_command, build_command]

    if sys.platform.startswith('win'):
        for vcvars_path in _get_vcvars_paths():
            if not vcvars_path.is_file():
                continue
            commands.insert(0, f'"{vcvars_path}" amd64')
            break
        else:
            raise Exit(
                'Error: Commands to configure MSVC environment variables not found.',
                code=1,
            )

    os.chdir(build_dir)
    c.run("&&".join(commands))
