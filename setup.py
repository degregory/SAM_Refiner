from distutils.core import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(  name='SAM_Refiner',
        version='1.0',
        py_modules=['SAM_Refiner'],
        description='A program for gathering variant information from a SAM formated files',
        long_description=long_description,
        long_description_content_type="text/markdown",
        author='Devon A. Gregory',
        author_email='gregoryde@missourie.edu',
        url='https://github.com/degregory/SAM_Refiner',
        license='GPL-3.0',
        entry_points={
            'console_scripts': [
                'SAM_Refiner = SAM_Refiner:main',
            ],
        },
      )