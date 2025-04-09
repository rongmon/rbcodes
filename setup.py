from setuptools import setup, find_packages

setup(
    name="rbcodes",
    version="0.1.0",
    description="A package for analyzing astronomical spectroscopic data",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/rbcodes",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
    "astropy>=5.3.3",
    "matplotlib==3.5.0",
    "numpy==1.22.3",
    "pandas==1.3.5",
    #"pyqt==5.15.7",
    "scipy==1.7.3",
    "scikit-learn>=1.5.0",
    ],
    python_requires=">=3.9",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)