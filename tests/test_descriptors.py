import numpy as np
import pyshot
import pytest

@pytest.mark.parametrize("radius", [0.1, 1., 2.])
@pytest.mark.parametrize("n_bins", [10, 20, 100])
@pytest.mark.parametrize("min_neighbors", [3, 5, 6])
@pytest.mark.parametrize("double_volumes_sectors", [True, False])
@pytest.mark.parametrize("use_interpolation", [True, False])
@pytest.mark.parametrize("use_normalization", [True, False])
def test_descriptors(
    radius,
    n_bins,
    min_neighbors,
    double_volumes_sectors,
    use_interpolation,
    use_normalization,
):

    verts = np.random.rand(500, 3)
    faces = np.random.choice(500, size=(1000, 3), replace=True)

    shot_features = pyshot.get_descriptors(
                verts,
                faces,
                radius=radius,
                local_rf_radius=radius,
                min_neighbors=min_neighbors,
                n_bins=n_bins,
                double_volumes_sectors=double_volumes_sectors,
                use_interpolation=use_interpolation,
                use_normalization=use_normalization,
    )

    # Coefficients must be finite
    assert np.isfinite(shot_features).all()

    n_features = 16 * (n_bins + 1) * (double_volumes_sectors + 1)
    assert shot_features.shape[1] == n_features