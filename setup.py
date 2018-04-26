from distutils.core import setup

setup(
    name='MESA2GADGET',
    packages=['MESA2GADGET', ],
    package_dir={'MESA2GADGET': 'src'}
    verstion='0.1.0',
    description='TBD',
    author='Meridith Joyce',
    author_email='TBD',
    url='https://github.com/mjoyceGR/MESA2GADGET',
    download_url='https://github.com/mjoyceGR/MESA2GADGET/archive/0.1.0.tar.gz',
    keywords=['mesa', 'gadget', 'astronomy'],
    classifiers=[],
    install_requires=[
        'h5py',
        'scipy',
        'healpy',
        'matplotlib'
    ]
)
