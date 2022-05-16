import sys, setuptools

if sys.version_info[0] < 3:
    sys.stdout.write("honto requires Python 3 or higher. Your Python your current Python version is {}.{}.{}"
                     .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))

setuptools.setup(name="honto",
                 version="0.1.0",
                 author="Fabio Cumbo",
                 author_email="fabio.cumbo@gmail.com",
                 url="http://github.com/fabio-cumbo/honto",
                 license="LICENSE",
                 packages=setuptools.find_packages(),
                 entry_points={
                     "console_scripts": ["honto = honto.honto:main"]
                 },
                 description="A novel method for assessing and measuring homophily in networks",
                 long_description=open("README.md").read(),
                 long_description_content_type="text/markdown",
                 install_requires=[
                     "tqdm",
                     "pandas",
                     "networkx",
                     "seaborn",
                     "matplotlib"
                 ],
                 zip_safe=False,
                 classifiers=[
                    "Programming Language :: Python :: 3",
                    "License :: OSI Approved :: MIT License",
                    "Operating System :: OS Independent",
                 ])