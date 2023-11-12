"""
Utility functions for

- building and importing modules on test time, using a temporary location
- detecting if compilers are present
- determining paths to tests

"""
import glob
import os
import sys
import subprocess
import tempfile
import shutil
import atexit
import textwrap
import re
import pytest
import contextlib
import numpy

from pathlib import Path
from numpy._utils import asunicode
from numpy.testing import temppath, IS_WASM
from importlib import import_module
from numpy.f2py._backends._meson import MesonBackend

#
# Maintaining a temporary module directory
#

_module_dir = None
_module_num = 5403

if sys.platform == "cygwin":
    NUMPY_INSTALL_ROOT = Path(__file__).parent.parent.parent
    _module_list = list(NUMPY_INSTALL_ROOT.glob("**/*.dll"))


def _cleanup():
    global _module_dir
    if _module_dir is not None:
        try:
            sys.path.remove(_module_dir)
        except ValueError:
            pass
        try:
            shutil.rmtree(_module_dir)
        except OSError:
            pass
        _module_dir = None


def get_module_dir():
    global _module_dir
    if _module_dir is None:
        _module_dir = tempfile.mkdtemp()
        atexit.register(_cleanup)
        if _module_dir not in sys.path:
            sys.path.insert(0, _module_dir)
    return _module_dir


def get_temp_module_name():
    # Assume single-threaded, and the module dir usable only by this thread
    global _module_num
    get_module_dir()
    name = "_test_ext_module_%d" % _module_num
    _module_num += 1
    if name in sys.modules:
        # this should not be possible, but check anyway
        raise RuntimeError("Temporary module name already in use.")
    return name


def _memoize(func):
    memo = {}

    def wrapper(*a, **kw):
        key = repr((a, kw))
        if key not in memo:
            try:
                memo[key] = func(*a, **kw)
            except Exception as e:
                memo[key] = e
                raise
        ret = memo[key]
        if isinstance(ret, Exception):
            raise ret
        return ret

    wrapper.__name__ = func.__name__
    return wrapper


#
# Building modules
#


@_memoize
def build_module(source_files, options=[], skip=[], only=[], module_name=None):
    """
    Compile and import a f2py module, built from the given files.

    """

    code = f"import sys; sys.path = {sys.path!r}; import numpy.f2py; numpy.f2py.main()"

    d = get_module_dir()

    # Copy files
    dst_sources = []
    f2py_sources = []
    for fn in source_files:
        if not os.path.isfile(fn):
            raise RuntimeError("%s is not a file" % fn)
        dst = os.path.join(d, os.path.basename(fn))
        shutil.copyfile(fn, dst)
        dst_sources.append(dst)

        base, ext = os.path.splitext(dst)
        if ext in (".f90", ".f", ".c", ".pyf"):
            f2py_sources.append(dst)

    assert f2py_sources

    # Prepare options
    if module_name is None:
        module_name = get_temp_module_name()
    f2py_opts = ["-c", "-m", module_name] + options + f2py_sources
    f2py_opts += ["--backend", "meson"]
    if skip:
        f2py_opts += ["skip:"] + skip
    if only:
        f2py_opts += ["only:"] + only

    # Build
    cwd = os.getcwd()
    try:
        os.chdir(d)
        cmd = [sys.executable, "-c", code] + f2py_opts
        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        out, err = p.communicate()
        if p.returncode != 0:
            raise RuntimeError("Running f2py failed: %s\n%s" %
                               (cmd[4:], asunicode(out)))
    finally:
        os.chdir(cwd)

        # Partial cleanup
        for fn in dst_sources:
            os.unlink(fn)

    # Rebase (Cygwin-only)
    if sys.platform == "cygwin":
        # If someone starts deleting modules after import, this will
        # need to change to record how big each module is, rather than
        # relying on rebase being able to find that from the files.
        _module_list.extend(
            glob.glob(os.path.join(d, "{:s}*".format(module_name)))
        )
        subprocess.check_call(
            ["/usr/bin/rebase", "--database", "--oblivious", "--verbose"]
            + _module_list
        )

    # Import
    return import_module(module_name)


@_memoize
def build_code(source_code,
               options=[],
               skip=[],
               only=[],
               suffix=None,
               module_name=None):
    """
    Compile and import Fortran code using f2py.

    """
    if suffix is None:
        suffix = ".f"
    with temppath(suffix=suffix) as path:
        with open(path, "w") as f:
            f.write(source_code)
        return build_module([path],
                            options=options,
                            skip=skip,
                            only=only,
                            module_name=module_name)


#
# Check if compilers are available at all...
#

def check_language(lang, code_snippet=None):
    tmpdir = tempfile.mkdtemp()
    try:
        meson_file = os.path.join(tmpdir, "meson.build")
        with open(meson_file, "w") as f:
            f.write("project('check_compilers')\n")
            f.write(f"add_languages('{lang}')\n")
            if code_snippet:
                f.write(f"{lang}_compiler = meson.get_compiler('{lang}')\n")
                f.write(f"{lang}_code = '''{code_snippet}'''\n")
                f.write(
                    f"_have_{lang}_feature ="
                    f"{lang}_compiler.compiles({lang}_code,"
                    f" name: '{lang} feature check')\n"
                )
        runmeson = subprocess.run(
            ["meson", "setup", "btmp"],
            check=False,
            cwd=tmpdir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if runmeson.returncode == 0:
            return True
        else:
            return False
    finally:
        shutil.rmtree(tmpdir)
    return False

fortran77_code = '''
C Example Fortran 77 code
      PROGRAM HELLO
      PRINT *, 'Hello, Fortran 77!'
      END
'''

fortran90_code = '''
! Example Fortran 90 code
program hello90
  type :: greeting
    character(len=20) :: text
  end type greeting

  type(greeting) :: greet
  greet%text = 'hello, fortran 90!'
  print *, greet%text
end program hello90
'''

# Dummy class for caching relevant checks
class CompilerChecker:
    def __init__(self):
        self.compilers_checked = False
        self.has_c = False
        self.has_f77 = False
        self.has_f90 = False

    def check_compilers(self):
        if not self.compilers_checked:
            self.has_c = check_language('c')
            self.has_f77 = check_language('fortran', fortran77_code)
            self.has_f90 = check_language('fortran', fortran90_code)
            self.compilers_checked = True

checker = CompilerChecker()
checker.check_compilers()

def has_c_compiler():
    return checker.has_c

def has_f77_compiler():
    return checker.has_f77

def has_f90_compiler():
    return checker.has_f90

#
# Building with distutils
#


@_memoize
def build_module_distutils(source_files, config_code, module_name, **kw):
    """
    Build a module via distutils and import it.

    """
    d = get_module_dir()

    # Copy files
    dst_sources = []
    for fn in source_files:
        if not os.path.isfile(fn):
            raise RuntimeError("%s is not a file" % fn)
        dst = os.path.join(d, os.path.basename(fn))
        shutil.copyfile(fn, dst)
        dst_sources.append(dst)

    # Build script
    config_code = textwrap.dedent(config_code).replace("\n", "\n    ")

    code = fr"""
import os
import sys
sys.path = {repr(sys.path)}

def configuration(parent_name='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('', parent_name, top_path)
    {config_code}
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
    """
    script = os.path.join(d, get_temp_module_name() + ".py")
    dst_sources.append(script)
    with open(script, "wb") as f:
        f.write(code.encode('latin1'))

    # Build
    cwd = os.getcwd()
    try:
        os.chdir(d)
        cmd = [sys.executable, script, "build_ext", "-i"]
        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)
        out, err = p.communicate()
        if p.returncode != 0:
            raise RuntimeError("Running distutils build failed: %s\n%s" %
                               (cmd[4:], asstr(out)))
    finally:
        os.chdir(cwd)

        # Partial cleanup
        for fn in dst_sources:
            os.unlink(fn)

    return import_module(module_name)


#
# Building with meson
#


class SimplifiedMesonBackend(MesonBackend):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def compile(self):
        self.write_meson_build(self.build_dir)
        self.run_meson(self.build_dir)


def build_meson(source_files, module_name, **kwargs):
    """
    Build a module via Meson and import it.
    """
    build_dir = tempfile.mkdtemp()

    # Initialize the MesonBackend
    backend = SimplifiedMesonBackend(
        modulename=module_name,
        sources=source_files,
        extra_objects=kwargs.get("extra_objects", []),
        build_dir=build_dir,
        include_dirs=kwargs.get("include_dirs", []),
        library_dirs=kwargs.get("library_dirs", []),
        libraries=kwargs.get("libraries", []),
        define_macros=kwargs.get("define_macros", []),
        undef_macros=kwargs.get("undef_macros", []),
        f2py_flags=kwargs.get("f2py_flags", []),
        sysinfo_flags=kwargs.get("sysinfo_flags", []),
        fc_flags=kwargs.get("fc_flags", []),
        flib_flags=kwargs.get("flib_flags", []),
        setup_flags=kwargs.get("setup_flags", []),
        remove_build_dir=kwargs.get("remove_build_dir", False),
        extra_dat=kwargs.get("extra_dat", {}),
    )

    # Compile the module
    backend.compile()

    # Import the compiled module
    sys.path.insert(0, f"{build_dir}/{backend.meson_build_dir}")
    return import_module(module_name)


#
# Unittest convenience
#


class F2PyTest:
    code = None
    sources = None
    options = []
    skip = []
    only = []
    suffix = ".f"
    module = None

    @property
    def module_name(self):
        cls = type(self)
        return f'_{cls.__module__.rsplit(".",1)[-1]}_{cls.__name__}_ext_module'

    def setup_method(self):
        if sys.platform == "win32":
            pytest.skip("Fails with MinGW64 Gfortran (Issue #9673)")

        if self.module is not None:
            return

        # Check compiler availability first
        if not has_c_compiler():
            pytest.skip("No C compiler available")

        codes = []
        if self.sources:
            codes.extend(self.sources)
        if self.code is not None:
            codes.append(self.suffix)

        needs_f77 = False
        needs_f90 = False
        needs_pyf = False
        for fn in codes:
            if str(fn).endswith(".f"):
                needs_f77 = True
            elif str(fn).endswith(".f90"):
                needs_f90 = True
            elif str(fn).endswith(".pyf"):
                needs_pyf = True
        if needs_f77 and not has_f77_compiler():
            pytest.skip("No Fortran 77 compiler available")
        if needs_f90 and not has_f90_compiler():
            pytest.skip("No Fortran 90 compiler available")
        if needs_pyf and not (has_f90_compiler() or has_f77_compiler()):
            pytest.skip("No Fortran compiler available")

        # Build the module
        if self.code is not None:
            self.module = build_code(
                self.code,
                options=self.options,
                skip=self.skip,
                only=self.only,
                suffix=self.suffix,
                module_name=self.module_name,
            )

        if self.sources is not None:
            self.module = build_module(
                self.sources,
                options=self.options,
                skip=self.skip,
                only=self.only,
                module_name=self.module_name,
            )


#
# Helper functions
#


def getpath(*a):
    # Package root
    d = Path(numpy.f2py.__file__).parent.resolve()
    return d.joinpath(*a)


@contextlib.contextmanager
def switchdir(path):
    curpath = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(curpath)
