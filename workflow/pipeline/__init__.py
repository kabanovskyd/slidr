import tomllib
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("slidr")
except PackageNotFoundError:
    # not pip/uv-installed with registered metadata (e.g. running uninstalled) -
    # fall back to reading the version directly out of pyproject.toml
    pyproject_path = Path(__file__).parent.parent.parent / "pyproject.toml"
    if pyproject_path.is_file():
        with open(pyproject_path, 'rb') as f:
            __version__ = tomllib.load(f).get("project", {}).get("version", "unknown")
    else:
        __version__ = "unknown"