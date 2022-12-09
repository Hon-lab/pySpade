from setuptools import setup, find_packages
from pySpade import __version__

with open(file='README.md', mode='r') as fh:
    long_description = fh.read()

    setup(
    name='pySpade',
    version=__version__,
    author='Yihan Wang',
    author_email='Yihan.Wang@UTSouthwestern.edu',
    description='Single cell Perturbation Analysis of Differential Expression',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/yihan1119/pySpade',
    packages=find_packages(exclude=('docs', 'tests')),
    package_dir={'pySpade': 'pySpade'},
    package_data={'pySpade': ['*.txt']},
    include_package_data=True,
    classifiers=[
                'Programming Language :: Python :: 3',
                'License :: OSI Approved :: MIT License',
                'Operating System :: OS Independent'
    ],
    python_requires='>=3.7',
    entry_points={
        'console_scripts': ['pySpade=pySpade.__main__:main'],
    }
)
