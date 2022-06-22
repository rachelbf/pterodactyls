# pterodactyls (Under Construction)
<a href="https://zenodo.org/badge/latestdoi/498841506"><img src="https://zenodo.org/badge/498841506.svg" alt="DOI"></a>

See [Fernandes R. B. et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv220603989F/abstract) for more information about ``pterodactyls``.

Prerequisites
-------------
1. Python = 3.7 (Might work with lower versions of Python but not higher due to tensorflow incompatibilty)
2. eleanor: https://github.com/afeinstein20/eleanor (edit max_sectors.py to 26 for primary mission data only)
3. wotan: https://github.com/hippke/wotan/blob/master/README.md
4. transitleastsquares: https://github.com/hippke/tls
5. triceratops: https://github.com/stevengiacalone/triceratops
6. EDIVetter Unplugged: https://github.com/jonzink/EDI_Vetter_unplugged
7. Vetting (for centroid test): https://pypi.org/project/vetting/
8. tabulate: https://pypi.org/project/tabulate/
9. lightkurve: https://docs.lightkurve.org


Usage
-------------
``pterodactyls`` can be easily used with jupyter notebook (with Python 3.7). See the example notebook for a brief tutorial.

## Attribution
If you use ``pterodactyls``, please cite [Fernandes R. B. et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv220603989F/abstract)
```
@ARTICLE{2022arXiv220603989F,
       author = {{Fernandes}, Rachel B. and {Mulders}, Gijs D. and {Pascucci}, Ilaria and {Bergsten}, Galen J. and {Koskinen}, Tommi T. and {Hardegree-Ullman}, Kevin K. and {Pearson}, Kyle A. and {Giacalone}, Steven and {Zink}, Jon and {Ciardi}, David R. and {O'Brien}, Patrick},
        title = "{Pterodactyls: A Tool to Uniformly Search and Vet for Young Transiting Planets In TESS Primary Mission Photometry}",
      journal = {arXiv e-prints},
     keywords = {Astrophysics - Earth and Planetary Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2022,
        month = jun,
          eid = {arXiv:2206.03989},
        pages = {arXiv:2206.03989},
archivePrefix = {arXiv},
       eprint = {2206.03989},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022arXiv220603989F},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
as well as the code:
```
@software{rachel_fernandes_2022_6667960,
  author       = {Rachel Fernandes},
  title        = {{pterodactyls: A Tool to Uniformly Search and Vet 
                   for Young Transiting Planets In TESS Primary
                   Mission Photometry}},
  month        = jun,
  year         = 2022,
  publisher    = {Zenodo},
  version      = {v0.1},
  doi          = {10.5281/zenodo.6667960},
  url          = {https://doi.org/10.5281/zenodo.6667960}
}
```
