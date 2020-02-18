
"""
MIT LICENSE

Copyright (c) 2019 @ZQWEI-Tech  Quix

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import setuptools

with open("README.MD", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rocketsimu",
    version="1.0.5",
    author="Quix Fan @ZQWEI-Tech",
    author_email="qxdnfsy@126.com",
    description="Simple tool for rocket simulation by python, the code were based on the original SRD  (Simple Rccket Designer)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://zqwei-quix.ukgsdn.co.uk/simple-rocket-designer/",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
        'scipy',
    ],
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)