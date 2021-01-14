import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="COT",
    version="1.0.0",
    author="Yingzhou Lu, Chiung-Ting Wu",
    author_email="lyz66@vt.edu, ctwu@vt.edu",
    description="Cosine based One-sample Test (COT) is an accurate and efficient method to detect subtype-specific marker genes (SMG) among many subtypes.",
    long_description=long_description,
    url="https://github.com/MintaYLu/COT",
    packages="COT",
    install_requires=[
       'numpy>=1.19.5',
       'pandas>=1.2.0',
       'scipy>=1.3.1',
       'statsmodels>=0.10.1',
       'scikit-learn>=0.21.3',
       'seaborn>=0.11.1',
       'matplotlib>=3.1.1'
    ],
    python_requires='>=3.7.4',
)