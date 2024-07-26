from setuptools import setup

def read_requirements():
    with open('requirements.txt') as req_file:
        return req_file.read().splitlines()

setup(
    name='graphpack',
    version='0.1.2',
    author='Bottazzi Daniele',
    author_email='daniele.bottazzi@mail.polimi.it',
    description='A Python tool to perform graph compression and visualization',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    #url='https://github.com/yourusername/graphpack',
    license='MIT',

    packages=['graphpack', 'graphpack.demo', 'graphpack.config'],
    package_dir={'': 'src'},

    install_requires=read_requirements(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',

    entry_points={
        'console_scripts': [
            'graphpack=graphpack.compression:main',
            'gp-demo=graphpack.demo.demo:main',
            'gp-sankey=graphpack.demo.sankey:main',
        ],
    },

    include_package_data=True,
    package_data={
        'graphpack': [
            'config/*',
            '../../data/input/*.txt'
        ],
    },
)
