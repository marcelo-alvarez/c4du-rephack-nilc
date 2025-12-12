"""Tests for ILC component separation module."""
import numpy as np
import pytest

# Note: pixell cannot be installed on Windows, so we'll mock it for testing
try:
    from pixell import enmap
    PIXELL_AVAILABLE = True
except ImportError:
    PIXELL_AVAILABLE = False
    # Create a mock enmap for testing when pixell is not available
    class MockEnmap:
        @staticmethod
        def zeros(shape, wcs=None):
            arr = np.zeros(shape)
            if hasattr(arr, '__array_finalize__'):
                arr.wcs = wcs
            return arr

        @staticmethod
        def zeros_like(template):
            return np.zeros_like(template)

        @staticmethod
        def ndmap(arr, wcs):
            if hasattr(arr, '__array_finalize__'):
                arr.wcs = wcs
            return arr

    enmap = MockEnmap()

from src.nilc.ilc.covariance import compute_inband_covariance
from src.nilc.ilc.weights import compute_ilc_weights, validate_weights
from src.nilc.ilc.bias_mitigation import (
    choose_bias_mitigation_strategy,
    apply_bias_mitigation,
)
from src.nilc.ilc.separation import apply_ilc_weights


def create_mock_map(shape=(100, 100), wcs=None):
    """Create a mock map for testing."""
    data = np.random.randn(*shape)
    if PIXELL_AVAILABLE:
        # Create proper CAR projection
        box = np.array([[-5, -5], [5, 5]]) * np.pi / 180  # 10x10 degree box
        return enmap.ndmap(data, enmap.geometry(shape, box)[1])
    else:
        # Return plain numpy array for testing without pixell
        return data


def create_mock_mask(shape=(100, 100), fill_fraction=0.8):
    """Create a mock mask."""
    mask = np.random.rand(*shape) < fill_fraction
    return mask.astype(float)


class TestCovariance:
    """Test covariance computation."""

    def test_compute_inband_covariance_shape(self):
        """Test that covariance has correct shape."""
        # Create two mock frequency maps
        map1 = create_mock_map()
        map2 = create_mock_map()
        needlet_maps = [map1, map2]
        mask = create_mock_mask()

        cov = compute_inband_covariance(needlet_maps, mask)

        assert cov.shape == (2, 2), "Covariance should be 2x2 for 2 frequencies"

    def test_covariance_symmetry(self):
        """Test that covariance matrix is symmetric."""
        map1 = create_mock_map()
        map2 = create_mock_map()
        needlet_maps = [map1, map2]
        mask = create_mock_mask()

        cov = compute_inband_covariance(needlet_maps, mask)

        assert np.allclose(cov, cov.T), "Covariance should be symmetric"

    def test_covariance_positive_diagonal(self):
        """Test that diagonal elements (variances) are positive."""
        map1 = create_mock_map()
        map2 = create_mock_map()
        needlet_maps = [map1, map2]
        mask = create_mock_mask()

        cov = compute_inband_covariance(needlet_maps, mask)

        assert np.all(np.diag(cov) > 0), "Variances should be positive"

    def test_covariance_empty_mask(self):
        """Test that empty mask raises error."""
        map1 = create_mock_map()
        map2 = create_mock_map()
        needlet_maps = [map1, map2]
        mask = np.zeros((100, 100))  # Empty mask

        with pytest.raises(ValueError, match="no valid pixels"):
            compute_inband_covariance(needlet_maps, mask)

    def test_covariance_with_perfect_correlation(self):
        """Test covariance when maps are identical."""
        # Create identical maps
        base_map = create_mock_map()
        needlet_maps = [base_map, base_map]
        mask = create_mock_mask()

        cov = compute_inband_covariance(needlet_maps, mask)

        # For identical maps, off-diagonal should equal diagonal
        assert np.isclose(cov[0, 1], cov[0, 0]), \
            "Identical maps should have equal diagonal and off-diagonal"


class TestWeights:
    """Test ILC weight computation."""

    def test_compute_ilc_weights_shape(self):
        """Test that weights have correct shape."""
        cov = np.array([[1.0, 0.3], [0.3, 1.0]])
        sed = np.array([1.0, 1.0])

        weights = compute_ilc_weights(cov, sed)

        assert weights.shape == (2,), "Weights should have length 2"

    def test_weights_unit_gain_cmb(self):
        """Test that CMB weights satisfy unit gain constraint."""
        cov = np.array([[1.0, 0.3], [0.3, 1.0]])
        cmb_sed = np.array([1.0, 1.0])

        weights = compute_ilc_weights(cov, cmb_sed)
        is_valid, gain = validate_weights(weights, cmb_sed)

        assert is_valid, f"Weights should satisfy unit gain, got {gain}"
        assert np.isclose(gain, 1.0), f"Gain should be 1.0, got {gain}"

    def test_weights_different_covariances(self):
        """Test that different covariances give different weights."""
        # Create two truly different covariances with different noise levels
        cov1 = np.array([[1.0, 0.1], [0.1, 2.0]])  # Different variances
        cov2 = np.array([[2.0, 0.3], [0.3, 1.0]])  # Different variances
        sed = np.array([1.0, 1.0])

        weights1 = compute_ilc_weights(cov1, sed)
        weights2 = compute_ilc_weights(cov2, sed)

        assert not np.allclose(weights1, weights2), \
            "Different covariances should yield different weights"

    def test_weights_with_singular_matrix(self):
        """Test handling of near-singular covariance."""
        # Create nearly singular matrix
        cov = np.array([[1.0, 0.9999], [0.9999, 1.0]])
        sed = np.array([1.0, 1.0])

        # Should not raise error due to regularization
        weights = compute_ilc_weights(cov, sed)
        is_valid, gain = validate_weights(weights, sed, tolerance=1e-3)

        assert is_valid, "Should handle near-singular matrices with regularization"

    def test_weights_dimension_mismatch(self):
        """Test error on dimension mismatch."""
        cov = np.array([[1.0, 0.3], [0.3, 1.0]])
        sed = np.array([1.0, 1.0, 1.0])  # Wrong dimension

        with pytest.raises(ValueError, match="must match"):
            compute_ilc_weights(cov, sed)


class TestBiasMitigation:
    """Test bias mitigation strategies."""

    def test_strategy_selection_large_scale(self):
        """Test that large scales use harmonic strategy."""
        ell_peak = 500
        strategy = choose_bias_mitigation_strategy(ell_peak, transition_ell=1000)

        assert strategy == 'harmonic', \
            f"ell={ell_peak} should use harmonic, got {strategy}"

    def test_strategy_selection_small_scale(self):
        """Test that small scales use donut strategy."""
        ell_peak = 2000
        strategy = choose_bias_mitigation_strategy(ell_peak, transition_ell=1000)

        assert strategy == 'donut', \
            f"ell={ell_peak} should use donut, got {strategy}"

    def test_strategy_at_transition(self):
        """Test strategy exactly at transition."""
        ell_peak = 1000
        strategy = choose_bias_mitigation_strategy(ell_peak, transition_ell=1000)

        assert strategy == 'donut', \
            f"ell={ell_peak} at transition should use donut, got {strategy}"

    @pytest.mark.skipif(not PIXELL_AVAILABLE, reason="Requires pixell")
    def test_apply_bias_mitigation_harmonic(self):
        """Test that harmonic bias mitigation returns filtered maps."""
        map1 = create_mock_map()
        map2 = create_mock_map()
        needlet_maps = [map1, map2]
        mask = create_mock_mask()
        ell_peak = 500  # Large scale -> harmonic

        processed, original, strategy = apply_bias_mitigation(
            needlet_maps, ell_peak, mask, transition_ell=1000
        )

        assert strategy == 'harmonic'
        assert len(processed) == 2
        assert processed is original  # For harmonic, they should be the same

    @pytest.mark.skipif(not PIXELL_AVAILABLE, reason="Requires pixell")
    def test_apply_bias_mitigation_donut(self):
        """Test that donut bias mitigation returns smoothed and original."""
        map1 = create_mock_map()
        map2 = create_mock_map()
        needlet_maps = [map1, map2]
        mask = create_mock_mask()
        ell_peak = 2000  # Small scale -> donut

        processed, original, strategy = apply_bias_mitigation(
            needlet_maps, ell_peak, mask, transition_ell=1000
        )

        assert strategy == 'donut'
        assert len(processed) == 2
        assert len(original) == 2
        # Processed should be smoothed versions
        assert processed is not original


class TestSeparation:
    """Test component separation."""

    def test_apply_ilc_weights_basic(self):
        """Test applying ILC weights to extract component."""
        # Create simple test maps
        map1 = np.ones((50, 50)) * 1.0
        map2 = np.ones((50, 50)) * 2.0
        frequency_maps = [map1, map2]
        weights = np.array([0.5, 0.5])

        result = apply_ilc_weights(frequency_maps, weights)

        expected = 0.5 * map1 + 0.5 * map2
        assert np.allclose(result, expected), "Weighted sum incorrect"

    def test_apply_ilc_weights_with_mask(self):
        """Test that mask is applied to output."""
        map1 = np.ones((50, 50))
        map2 = np.ones((50, 50)) * 2.0
        frequency_maps = [map1, map2]
        weights = np.array([0.5, 0.5])
        mask = np.zeros((50, 50))
        mask[:25, :] = 1.0  # Half masked

        result = apply_ilc_weights(frequency_maps, weights, mask)

        # Masked region should be zero
        assert np.all(result[25:, :] == 0), "Masked region should be zero"
        # Unmasked region should have values
        assert np.all(result[:25, :] != 0), "Unmasked region should have values"

    def test_apply_ilc_weights_dimension_mismatch(self):
        """Test error when weights don't match maps."""
        map1 = np.ones((50, 50))
        map2 = np.ones((50, 50))
        frequency_maps = [map1, map2]
        weights = np.array([0.5, 0.5, 0.0])  # Wrong dimension

        with pytest.raises(ValueError, match="must match"):
            apply_ilc_weights(frequency_maps, weights)

    def test_apply_ilc_weights_preserves_cmb(self):
        """Test that CMB signal is preserved with correct weights."""
        # Create maps with CMB signal
        cmb_signal = np.random.randn(50, 50)
        noise1 = np.random.randn(50, 50) * 0.1
        noise2 = np.random.randn(50, 50) * 0.1

        # Both frequencies see same CMB (achromatic)
        map1 = cmb_signal + noise1
        map2 = cmb_signal + noise2
        frequency_maps = [map1, map2]

        # Equal weights for achromatic component
        weights = np.array([0.5, 0.5])

        result = apply_ilc_weights(frequency_maps, weights)

        # Result should be close to CMB (averaging reduces noise)
        correlation = np.corrcoef(result.flatten(), cmb_signal.flatten())[0, 1]
        assert correlation > 0.9, "Should preserve CMB signal"


class TestIntegration:
    """Integration tests for full pipeline."""

    @pytest.mark.skipif(not PIXELL_AVAILABLE, reason="Requires pixell")
    def test_end_to_end_separation(self):
        """Test end-to-end component separation with mock data."""
        from src.nilc.ilc import separate_components

        # Create mock needlet-decomposed maps
        ell_peaks = [100, 200]
        needlet_maps = {}

        for ell_peak in ell_peaks:
            map1 = create_mock_map()
            map2 = create_mock_map()
            needlet_maps[ell_peak] = [map1, map2]

        mask = create_mock_mask()

        # Run separation
        separated = separate_components(
            needlet_maps, mask, ell_peaks, components=['cmb', 'ksz']
        )

        # Check output structure
        assert 'cmb' in separated
        assert 'ksz' in separated
        assert 'weights' in separated
        assert 'metadata' in separated

        # Check that each component has maps for each scale
        for ell_peak in ell_peaks:
            assert ell_peak in separated['cmb']
            assert ell_peak in separated['ksz']
            assert ell_peak in separated['weights']['cmb']
            assert ell_peak in separated['weights']['ksz']

        # Check metadata
        assert 'strategies' in separated['metadata']
        for ell_peak in ell_peaks:
            assert ell_peak in separated['metadata']['strategies']

    def test_separation_with_synthetic_components(self):
        """Test separation with known synthetic CMB and kSZ."""
        # Create synthetic components
        cmb = np.random.randn(50, 50)
        ksz = np.random.randn(50, 50) * 0.5
        noise1 = np.random.randn(50, 50) * 0.2
        noise2 = np.random.randn(50, 50) * 0.2

        # Both CMB and kSZ are achromatic
        map1 = cmb + ksz + noise1
        map2 = cmb + ksz + noise2

        # Simple ILC with equal weights should preserve total signal
        frequency_maps = [map1, map2]
        weights = np.array([0.5, 0.5])

        result = apply_ilc_weights(frequency_maps, weights)

        # Result should correlate with CMB+kSZ
        true_signal = cmb + ksz
        correlation = np.corrcoef(result.flatten(), true_signal.flatten())[0, 1]

        assert correlation > 0.8, \
            f"Should preserve achromatic signal, correlation={correlation}"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
