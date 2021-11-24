from setuptools import setup

setup(
        name='HippsDimes',
        version='1.0',
        py_modules=['HippsDimes'],
        install_requires=[
            'Click',
            'Numpy',
            'Scipy',
            'Pandas',
            'Tqdm',
            'Cooler',
            'rich',
            ],
        entry_points='''
            [console_scripts]
            HippsDimes=HippsDimes:main
            ''',
            )

