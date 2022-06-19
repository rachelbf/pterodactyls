# pterodactyls (Under Construction.....)

Python Tool for Exoplanets: Really Outstanding Detection and Assessment of Close-in Transits around Young Local Stars

See <a href="https://ui.adsabs.harvard.edu/abs/2022arXiv220603989F/abstract">Fernandes et al. 2022</a> for more information about this tool.

Usage
-------------

``pterodactyls`` can be easily used with jupyter notebook (with Python 3.7). See the notebook in the examples/ directory for a brief tutorial.

Prerequisites
-------------

``pterodactyls`` requires the following in order to work smoothly:
1. Python 3.7 (Can potentially work with versions lower than 3.7 but not higher due to Tensorflow compatibility issues)
2. eleanor: https://github.com/afeinstein20/eleanor (edit max_sectors.py to 26 for primary mission data only)
3. wotan: https://github.com/hippke/wotan/blob/master/README.md
4. transitleastsquares: https://github.com/hippke/tls
5. triceratops: https://github.com/stevengiacalone/triceratops
6. vetting (for centroid test): https://pypi.org/project/vetting/
7. EDIVetter Unplugged: https://github.com/jonzink/EDI_Vetter_unplugged
8. tabulate: https://pypi.org/project/tabulate/
9. lightkurve: https://docs.lightkurve.org

Attribution
-------------
If you use ``pterodactyls``, please cite both the paper and the code.

Paper citation::

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
