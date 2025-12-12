"""Pytest configuration and shared fixtures for NILC pipeline tests."""

import pytest
import numpy as np


@pytest.fixture
def sample_map_shape():
    """Standard test map shape (ny, nx)."""
    return (100, 100)


@pytest.fixture
def sample_frequencies():
    """ACT frequency channels in GHz."""
    return [90, 150]


@pytest.fixture
def sample_needlet_ells():
    """Needlet scale peaks for testing (subset of full pipeline)."""
    return [100, 200, 400, 800, 1600]
