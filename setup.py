from setuptools import setup, find_packages

setup(
    name='spinncer',
    version='0.0.1',
    packages=find_packages(),
    url='https://github.com/pabogdan/spinncer',
    license="GNU GPLv3.0",
    author='Petrut Antoniu Bogdan',
    author_email='petrut.bogdan@manchester.ac.uk',
    description='Simulating a large scale Cerebellum model on SpiNNaker',
    # Requirements
    dependency_links=[],

    install_requires=["numpy",
                      "scipy",
                      "brian2",
                      "spynnaker>= 1!5.0.0, < 1!6.0.0",
                      "pynn",
                      "matplotlib"],
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
