"""
Tests for ILC component separation module.

Tests cover covariance computation, weight calculation, bias mitigation,
and full component separation pipeline.
"""
import numpy as np
import pytest
from pixell import enmap
from src.nilc.ilc import (
    compute_covariance_map,
    compute_ilc_weights,
    compute_constrained_ilc_weights,
    get_cmb_response,
    get_ksz_response,
    get_tsz_response,
    apply_ilc_weights,
    separate_components_single_scale,
    invert_covariance,
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


class TestCovarianceComputation:
    """Test covariance matrix computation."""

    def test_covariance_shape(self):
        """Covariance map should have correct shape."""
        # Create mock needlet coefficients
        n_freq, n_y, n_x = 2, 10, 10
        needlet_coeffs = np.random.randn(n_freq, n_y, n_x)

        cov_map = compute_covariance_map(needlet_coeffs)

        assert cov_map.shape == (n_y, n_x, n_freq, n_freq)

    def test_covariance_symmetry(self):
        """Covariance matrix should be symmetric."""
        n_freq, n_y, n_x = 2, 10, 10
        needlet_coeffs = np.random.randn(n_freq, n_y, n_x)

        cov_map = compute_covariance_map(needlet_coeffs)

        # Check symmetry at each pixel
        for iy in range(n_y):
            for ix in range(n_x):
                C = cov_map[iy, ix, :, :]
                assert np.allclose(C, C.T)

    def test_covariance_positive_semidefinite(self):
        """Covariance matrix should be positive semi-definite."""
        n_freq, n_y, n_x = 2, 10, 10
        needlet_coeffs = np.random.randn(n_freq, n_y, n_x)

        cov_map = compute_covariance_map(needlet_coeffs)

        # Check positive semi-definite at a few random pixels
        for _ in range(5):
            iy, ix = np.random.randint(0, n_y), np.random.randint(0, n_x)
            C = cov_map[iy, ix, :, :]
            eigenvalues = np.linalg.eigvalsh(C)
            assert np.all(eigenvalues >= -1e-10)  # Allow small numerical error

    def test_covariance_inversion(self):
        """Test covariance matrix inversion."""
        n_freq, n_y, n_x = 2, 5, 5
        # Create well-conditioned covariance
        needlet_coeffs = np.random.randn(n_freq, n_y, n_x)
        cov_map = compute_covariance_map(needlet_coeffs)

        # Add regularization to diagonal
        for i in range(n_freq):
            cov_map[:, :, i, i] += 1.0

        cov_inv = invert_covariance(cov_map, regularization=1e-6)

        # Check that C @ C^(-1) ≈ I at a pixel
        iy, ix = 2, 2
        C = cov_map[iy, ix, :, :]
        C_inv = cov_inv[iy, ix, :, :]
        identity = C @ C_inv
        assert np.allclose(identity, np.eye(n_freq), atol=1e-5)


class TestILCWeights:
    """Test ILC weight computation."""

    def test_ilc_weights_shape(self):
        """ILC weights should have correct shape."""
        n_freq, n_y, n_x = 2, 10, 10

        # Create mock inverse covariance
        cov_inv = np.random.randn(n_y, n_x, n_freq, n_freq)

        # Component vector (CMB: achromatic)
        component_vector = np.ones(n_freq)

        weights = compute_ilc_weights(cov_inv, component_vector)

        assert weights.shape == (n_y, n_x, n_freq)

    def test_ilc_weights_normalization(self):
        """ILC weights should preserve the component (sum to 1 for achromatic)."""
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

    def test_constrained_ilc_shape(self):
        """Constrained ILC should have correct shape."""
        n_freq, n_y, n_x = 2, 5, 5
        n_comp = 2

        # Create mock inverse covariance
        cov_inv = np.zeros((n_y, n_x, n_freq, n_freq))
        for iy in range(n_y):
            for ix in range(n_x):
                cov_inv[iy, ix, :, :] = np.eye(n_freq)

        # Response matrix (CMB and kSZ, both achromatic)
        response_matrix = np.ones((n_freq, n_comp))

        weights = compute_constrained_ilc_weights(cov_inv, response_matrix)

        assert weights.shape == (n_y, n_x, n_freq, n_comp)


class TestComponentSeparation:
    """Test component separation functions."""

    def test_apply_weights(self):
        """Test applying ILC weights to data."""
        n_freq, n_y, n_x = 2, 10, 10

        # Create mock data (constant signal + noise)
        signal = 10.0
        needlet_coeffs = np.ones((n_freq, n_y, n_x)) * signal
        needlet_coeffs += np.random.randn(n_freq, n_y, n_x) * 0.1  # Add small noise

        # Create equal weights (should average the frequencies)
        ilc_weights = np.ones((n_y, n_x, n_freq)) / n_freq

        result = apply_ilc_weights(needlet_coeffs, ilc_weights)

        # Result should be close to the signal (average of frequencies)
        assert result.shape == (n_y, n_x)
        assert np.allclose(result, signal, atol=0.5)

    def test_single_scale_separation(self):
        """Test component separation for single needlet scale."""
        # Create simple test case
        n_freq, n_y, n_x = 2, 8, 8
        n_comp = 1

        # Create mock data (CMB signal)
        cmb_signal = np.random.randn(n_y, n_x)
        needlet_coeffs = np.zeros((n_freq, n_y, n_x))
        needlet_coeffs[0, :, :] = cmb_signal  # Same signal at both frequencies
        needlet_coeffs[1, :, :] = cmb_signal
        # Add independent noise
        needlet_coeffs += np.random.randn(n_freq, n_y, n_x) * 0.1

        # Create WCS (simple CAR projection)
        shape = (n_y, n_x)
        wcs = enmap.geometry(pos=(0, 0), shape=shape, res=1.0*np.pi/180)[1]

        # CMB response
        response_matrix = np.ones((n_freq, n_comp))

        # Separate (without bias mitigation for simplicity)
        components, weights_out = separate_components_single_scale(
            needlet_coeffs,
            ell_peak=1000,
            wcs=wcs,
            component_responses=response_matrix,
            mask=None,
            use_bias_mitigation=False
        )

        assert components.shape == (n_comp, n_y, n_x)
        assert weights_out.shape == (n_y, n_x, n_freq, n_comp)

        # Recovered CMB should be roughly correlated with input
        # (exact recovery not expected due to noise)
        recovered_cmb = components[0, :, :]
        correlation = np.corrcoef(cmb_signal.flatten(), recovered_cmb.flatten())[0, 1]
        assert correlation > 0.5  # Should have positive correlation


class TestIntegration:
    """Integration tests for full workflow."""

    def test_two_frequency_separation(self):
        """Test separating CMB from two ACT-like frequencies."""
        # Simulate simple CMB + noise at two frequencies
        n_y, n_x = 16, 16
        shape = (n_y, n_x)
        wcs = enmap.geometry(pos=(0, 0), shape=shape, res=1.0*np.pi/180)[1]

        # Create CMB signal
        cmb_map = np.random.randn(n_y, n_x)

        # Create "observations" at 90 and 150 GHz
        # CMB is achromatic, so same signal at both
        freq90 = cmb_map + np.random.randn(n_y, n_x) * 1.0  # Higher noise
        freq150 = cmb_map + np.random.randn(n_y, n_x) * 0.5  # Lower noise

        needlet_coeffs = np.array([freq90, freq150])

        # Response matrix for CMB only
        response_matrix = np.ones((2, 1))

        # Separate
        components, _ = separate_components_single_scale(
            needlet_coeffs,
            ell_peak=1000,
            wcs=wcs,
            component_responses=response_matrix,
            use_bias_mitigation=False
        )

        recovered_cmb = components[0]

        # ILC should do better than either individual frequency
        # (variance reduction through optimal weighting)
        ilc_noise = np.std(recovered_cmb - cmb_map)
        freq90_noise = np.std(freq90 - cmb_map)
        freq150_noise = np.std(freq150 - cmb_map)

        assert ilc_noise < freq90_noise  # Should beat noisier channel
        assert ilc_noise < freq150_noise  # Should beat even the better channel


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
