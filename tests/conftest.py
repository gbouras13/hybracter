"""
holds command line options for pytest

"""


def pytest_addoption(parser):
    parser.addoption("--run_mac", action="store_true")
