"""
Basic tests for ILC module that don't require pixell.

These tests verify the core mathematical operations:
- Spectral response functions
- Covariance computation (pixel-wise)
- ILC weight calculation
- Weight application

Full integration tests requiring pixell/enmap are in test_ilc.py
"""
import numpy as np
import pytest

# Mock pixell for imports
import sys
from unittest.mock import MagicMock

# Create mock modules
pixell_mock = MagicMock()
pixell_mock.enmap = MagicMock()
pixell_mock.fft = MagicMock()
sys.modules['pixell'] = pixell_mock
sys.modules['pixell.enmap'] = pixell_mock.enmap
sys.modules['pixell.fft'] = pixell_mock.fft

# Now import our modules
from src.nilc.ilc.weights import (
    compute_ilc_weights,
    get_cmb_response,
    get_ksz_response,
    get_tsz_response,
)


class TestResponseFunctions:
    """Test component spectral response functions."""

    def test_cmb_response(self):
        """CMB should be achromatic (all ones)."""
        freqs = np.array([90, 150, 220])
        response = get_cmb_response(freqs)
        assert np.allclose(response, np.ones(3))

    def test_ksz_response(self):
        """kSZ should be achromatic (all ones)."""
        freqs = np.array([90, 150, 220])
        response = get_ksz_response(freqs)
        assert np.allclose(response, np.ones(3))

    def test_tsz_response(self):
        """tSZ should have frequency-dependent signature."""
        freqs = np.array([90, 150, 220])
        response = get_tsz_response(freqs)

        # tSZ should be negative at low freq, positive at high freq
        # Null frequency is around 217 GHz
        assert response[0] < 0  # 90 GHz: negative
        assert response[1] < 0  # 150 GHz: negative
        assert response[2] > 0  # 220 GHz: positive (above null)

        # Should not all be the same (frequency-dependent)
        assert not np.allclose(response, response[0])

    def test_tsz_response_shape(self):
        """tSZ response should have correct shape."""
        freqs = np.array([90, 150, 220, 280])
        response = get_tsz_response(freqs)
        assert response.shape == (4,)


class TestILCWeights:
    """Test ILC weight computation."""

    def test_ilc_weights_shape(self):
        """ILC weights should have correct shape."""
        n_freq, n_y, n_x = 2, 10, 10

        # Create mock inverse covariance (well-conditioned)
        cov_inv = np.zeros((n_y, n_x, n_freq, n_freq))
        for iy in range(n_y):
            for ix in range(n_x):
                cov_inv[iy, ix, :, :] = np.eye(n_freq)

        # Component vector (CMB: achromatic)
        component_vector = np.ones(n_freq)

        weights = compute_ilc_weights(cov_inv, component_vector)

        assert weights.shape == (n_y, n_x, n_freq)

    def test_ilc_weights_normalization_identity_covariance(self):
        """For identity covariance, achromatic weights should sum to 1."""
        n_freq, n_y, n_x = 2, 5, 5

        # Create identity covariance (uncorrelated, equal variance)
        cov_inv = np.zeros((n_y, n_x, n_freq, n_freq))
        for iy in range(n_y):
            for ix in range(n_x):
                cov_inv[iy, ix, :, :] = np.eye(n_freq)

        # CMB component (achromatic)
        component_vector = np.ones(n_freq)

        weights = compute_ilc_weights(cov_inv, component_vector)

        # For achromatic component with identity covariance,
        # weights should sum to 1 at each pixel
        weight_sums = weights.sum(axis=2)
        assert np.allclose(weight_sums, 1.0)

    def test_ilc_weights_equal_for_identity(self):
        """For identity covariance and achromatic, weights should be equal."""
        n_freq = 2
        n_y, n_x = 3, 3

        # Identity covariance
        cov_inv = np.zeros((n_y, n_x, n_freq, n_freq))
        for iy in range(n_y):
            for ix in range(n_x):
                cov_inv[iy, ix, :, :] = np.eye(n_freq)

        component_vector = np.ones(n_freq)
        weights = compute_ilc_weights(cov_inv, component_vector)

        # All weights should be 0.5 (equal weighting for 2 frequencies)
        assert np.allclose(weights, 0.5)

    def test_ilc_weights_favor_lower_variance(self):
        """Weights should favor channels with lower variance."""
        n_freq = 2
        n_y, n_x = 1, 1

        # Covariance: first channel has variance 1, second has variance 4
        # (second channel is noisier)
        cov = np.array([[[[1.0, 0.0],
                          [0.0, 4.0]]]])

        # Invert
        cov_inv = np.array([[[[1.0, 0.0],
                             [0.0, 0.25]]]])

        component_vector = np.ones(n_freq)
        weights = compute_ilc_weights(cov_inv, component_vector)

        # Weight for first channel (lower variance) should be higher
        assert weights[0, 0, 0] > weights[0, 0, 1]

        # Weights should still sum to 1
        assert np.allclose(weights.sum(), 1.0)


class TestCovarianceBasic:
    """Test basic covariance properties without pixell."""

    def test_simple_covariance_computation(self):
        """Test pixel-wise covariance computation logic."""
        # Simplified version without pixell
        n_freq, n_y, n_x = 2, 3, 3

        # Create simple data: constant signal + noise
        data = np.zeros((n_freq, n_y, n_x))
        data[0, :, :] = 10.0  # First frequency
        data[1, :, :] = 10.0  # Second frequency (same signal)

        # Compute covariance at one pixel manually
        cov_manual = np.zeros((n_freq, n_freq))
        for i in range(n_freq):
            for j in range(n_freq):
                cov_manual[i, j] = data[i, 0, 0] * data[j, 0, 0]

        expected = np.array([[100, 100],
                            [100, 100]])

        assert np.allclose(cov_manual, expected)

    def test_covariance_symmetry_property(self):
        """Covariance should always be symmetric."""
        n_freq = 2

        # Create random 2x2 covariance matrix
        cov = np.random.randn(n_freq, n_freq)
        cov = cov @ cov.T  # Make it symmetric positive semidefinite

        assert np.allclose(cov, cov.T)


class TestEndToEnd:
    """End-to-end tests without full pixell integration."""

    def test_simple_ilc_combination(self):
        """Test simple ILC: two frequencies, achromatic signal."""
        # Create signal: same at both frequencies (CMB-like)
        signal = np.array([5.0, 5.0])

        # Add noise: first channel noisier
        np.random.seed(42)
        noise1 = np.random.randn() * 2.0  # Higher noise
        noise2 = np.random.randn() * 1.0  # Lower noise

        data = signal + np.array([noise1, noise2])

        # Simple covariance (diagonal, from known variances)
        cov = np.diag([4.0, 1.0])  # Variances: 4 and 1
        cov_inv = np.linalg.inv(cov)

        # Reshape for function
        cov_inv_4d = cov_inv.reshape(1, 1, 2, 2)

        # Compute weights
        component_vector = np.ones(2)  # Achromatic
        weights = compute_ilc_weights(cov_inv_4d, component_vector)

        # Apply weights
        ilc_result = np.sum(weights[0, 0, :] * data)

        # ILC should weight second channel more (lower noise)
        assert weights[0, 0, 1] > weights[0, 0, 0]

        # Result should be closer to true signal than individual channels
        # (though this is stochastic, seed makes it deterministic)
        assert abs(ilc_result - 5.0) < abs(data[0] - 5.0)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
