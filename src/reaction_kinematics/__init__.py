from importlib.metadata import PackageNotFoundError, version

from reaction_kinematics.reaction import Reaction

try:
    __version__ = version("reaction-kinematics")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = ["Reaction", "__version__"]
