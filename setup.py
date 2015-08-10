from ast import literal_eval
from distutils.core import setup


def get_version(source='src/phasing/__init__.py'):
    with open(source) as f:
        for line in f:
            if line.startswith('__version__'):
                return literal_eval(line.split('=')[-1].lstrip())
    raise ValueError("__version__ not found")


setup(
    name='phasing',
    version=get_version(),
    author='Nicholas Harding',
    author_email='njh@well.ox.ac.uk',
    package_dir={'': 'src'},
    packages=['phasing'],
    url='https://github.com/hardingnj/phasing',
    license='MIT License',
    description='General code to assist in phasing analysis',
    long_description=open('README.md').read(),
    classifiers=['Intended Audience :: Developers',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python',
                 'Topic :: Software Development :: Libraries :: Python Modules'
                 ],
    scripts=['bin/hdf5_2_vcf.py',
             'bin/select_truth_variants.py', 'bin/report_ME.py',
             'bin/identify_roh.py',
             'bin/evaluate_phasing_perf.py'])
