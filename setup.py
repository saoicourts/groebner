import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="groebner-nicocourts",
    version="0.0.1",
    author="Nico Courts",
    author_email="nico@nicocourts.com",
    description="A package implementing polynomials and Groebner bases",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NicoCourts/groebner",
    project_urls={
        "Bug Tracker": "https://github.com/NicoCourts/groebner/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
)