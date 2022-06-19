from setuptools import setup, find_packages
from pterodactyls_30 import __version__

setup(
    name='pterodactyls',
    version=__version__,
    description='Python Tool for Exoplanets: Really Outstanding Detection and Assessment of Close-in Transits around Young Local Stars',
    url='https://github.com/rachelbf/pterodactyls',
    author='Rachel B. Fernandes',
    author_email='rachelbf@lpl.arizona.edu',
    license='MIT',
    zip_safe=False,
    packages = find_packages(),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Astronomy'
        ],
    install_requires=[
        'numpy', 'matplotlib', 'pandas' ,'eleanor', 'wotan', 'transitleastsquares', 'triceratops', 'vetting'
        ]
)