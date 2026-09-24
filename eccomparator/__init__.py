import warnings

warnings.filterwarnings(
    "ignore",
    message="pkg_resources is deprecated",
    category=UserWarning,
)

__version__ = "0.1.3"
