from setuptools import setup, find_packages


install_requires = [
    "numpy > 1.13, < 1.20; python_version == '3.6'",
    "numpy > 1.13, < 1.21; python_version == '3.7'",
    "numpy; python_version >= '3.8'",
    "scipy >= 0.16.0, < 1.6; python_version == '3.6'",
    "scipy >= 0.16.0, < 1.8; python_version == '3.7'",
    "scipy >= 0.16.0; python_version >= '3.8'",
    "brian2",
    "neo==0.7.1",
    "pynn==0.9.5",
    "h5py",
    "matplotlib",
    "argparse",
    "elephant==0.6.2",
    "colorama",
    "pandas==0.23.4",
    "hilbertcurve"
]

setup(
    name='spinncer',
    version='1.0.0',
    packages=find_packages(),
    url='https://github.com/pabogdan/spinncer',
    license="GNU GPLv3.0",
    author='Petrut Antoniu Bogdan',
    author_email='petrut.bogdan@manchester.ac.uk',
    description='Simulating a large scale Cerebellum model on SpiNNaker',
    # Requirements
    install_requires=install_requires,
    classifiers=[
        "Development Status :: 3 - Alpha",

        "Intended Audience :: Science/Research",

        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",

        "Programming Language :: Python :: 3"
        "Programming Language :: Python :: 3.7"

        "Topic :: Scientific/Engineering",
    ]
)
