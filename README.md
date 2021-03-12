# `pyshot`: a Python interface for the SHOT descriptor

The original implementation originates from the author of [1], [2], [3].

It has been marginally modified to interface the computation of the SHOT
descriptors with other languages.

This project makes use of Cython to bind the implementation to a simple Python
interface.

## Installation step

 - Clone the repository
 - Install the project using the source
```bash
git clone <this repo url>
cd pyshot
# Manage your virtual env here.
# You can use, but don't need a conda env.
pip install .
```

### ⚠ C and C++ dependencies 
You need `Eigen3`, `FLANN` and `LZ4` to be able to compile the C++
implementation. 

Please refer to your distribution's package manager.

### Usage

Typical snippet:

```python
import pyshot
import numpy as np

### ...

vertices: np.array =  # ... a np.array of shape (n, 3)
faces: np.array =  # ... a np.array of shape (m, 3)

# a np.array of shape (n, n_descr)
shot_descrs: np.array = pyshot.get_descriptors(
	vertices,
	faces,
	radius=100,
	local_rf_radius=100,
	# The following parameters are optional
	min_neighbors=3,
	n_bins=20,
	double_volumes_sectors=True,
	use_interpolation=True,
	use_normalization=True,
)
```

See [`example_pyshot.py`](./example_pyshot.py).


## Installation step for development mode

 - Clone the repository.
 - Install the project using the source in [editable mode](https://packaging.python.org/guides/distributing-packages-using-setuptools/#working-in-development-mode)
```bash
git clone <this repo url>
cd pyshot
# Manage your virtual env here.
# You can use, but don't need a conda env.
pip install --editable . -v
```

### References

If you use this implementation in your work, please cite the following publications:

> [1] 	F. Tombari *, S. Salti *, L. Di Stefano, "Unique Signatures of Histograms for Local Surface Description",11th European Conference on Computer Vision (ECCV), September 5-11, Hersonissos, Greece, 2010.

> [2] 	F. Tombari, S. Salti, L. Di Stefano, "A combined texture-shape descriptor for enhanced 3D feature matching", IEEE International Conference on Image Processing (ICIP), September 11-14, Brussels, Belgium, 2011. [PDF] 

> [3] 	S. Salti, F. Tombari, L. Di Stefano, "SHOT: Unique Signatures of Histograms for Surface and Texture Description", Computer Vision and Image Understanding, May, 2014. [PDF] 
  	(* indicates equal contribution)

using:

```bibtex
@inproceedings{10.5555/1927006.1927035,
    author  = {Federico Tombari and
               Samuele Salti and
               Luigi di Stefano},
    title = {Unique Signatures of Histograms for Local Surface Description},
    year = {2010},
    isbn = {364215557X},
    publisher = {Springer-Verlag},
    address = {Berlin, Heidelberg},
    abstract = {This paper deals with local 3D descriptors for surface matching. First, we categorize existing methods into two classes: Signatures and Histograms. Then, by discussion and experiments alike, we point out the key issues of uniqueness and repeatability of the local reference frame. Based on these observations, we formulate a novel comprehensive proposal for surface representation, which encompasses a new unique and repeatable local reference frame as well as a new 3D descriptor. The latter lays at the intersection between Signatures and Histograms, so as to possibly achieve a better balance between descriptiveness and robustness. Experiments on publicly available datasets as well as on range scans obtained with Spacetime Stereo provide a thorough validation of our proposal.},
    booktitle = {Proceedings of the 11th European Conference on Computer Vision Conference on Computer Vision: Part III},
    pages = {356–369},
    numpages = {14},
    location = {Heraklion, Crete, Greece},
    series = {ECCV'10}
}

@inproceedings{DBLP:conf/icip/TombariSS11,
    author    = {Federico Tombari and
               Samuele Salti and
               Luigi di Stefano},
    editor    = {Beno{\^{\i}}t Macq and
               Peter Schelkens},
    title     = {A combined texture-shape descriptor for enhanced 3D feature matching},
    booktitle = {18th {IEEE} International Conference on Image Processing, {ICIP} 2011,
               Brussels, Belgium, September 11-14, 2011},
    pages     = {809--812},
    publisher = {{IEEE}},
    year      = {2011},
    url       = {https://doi.org/10.1109/ICIP.2011.6116679},
    doi       = {10.1109/ICIP.2011.6116679},
    timestamp = {Wed, 16 Oct 2019 14:14:52 +0200},
    biburl    = {https://dblp.org/rec/conf/icip/TombariSS11.bib},
    bibsource = {dblp computer science bibliography, https://dblp.org}
}

@article{SALTI2014251,
    author    = {Federico Tombari and
               Samuele Salti and
               Luigi di Stefano},
    title = {SHOT: Unique signatures of histograms for surface and texture description},
    journal = {Computer Vision and Image Understanding},
    volume = {125},
    pages = {251-264},
    year = {2014},
    issn = {1077-3142},
    doi = {https://doi.org/10.1016/j.cviu.2014.04.011},
    url = {https://www.sciencedirect.com/science/article/pii/S1077314214000988},
    keywords = {Surface matching, 3D descriptors, Object recognition, 3D reconstruction},
}
```

## License ([`fedassa/SHOT`](https://github.com/fedassa/SHOT) for reference)

SHOT is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

SHOT is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with SHOT. If not, see http://www.gnu.org/licenses/.