---
layout: default
title: Local testing
excerpt:
nav_order: 2
---

# Testing the documentation locally

Before commiting any changes to github we should test them locally.
We'll need ruby and many required gems (including jekyll and github-pages) installed in a virtual environment using conda.

```
conda create -n jekyll-test
conda activate jekyll-test

# apparently ruby packages at conda-forge are messed, this should work
#     conda install ruby rb-bundler
# but in my case I needed
conda install ruby rb-bundler c-compiler compilers cxx-compiler

# now change to the directory with documentation: `<somedir>/TDFstat/docs` and run
bundle install
# (this will install all the gems)
```

Having the virtual environment ready it is very easy to view the documentation locally.
Just change to `<somedir>/TDFstat/docs` directory and run:
```
bundle exec jekyll serve -l

```
This command generates static pages with jekyll and starts local web server hosting them at address listed in the output, usually http://127.0.0.1:4000 . Open this location in any web browser. \\
To facilitate editing, if any file in the docs directory is modified the pages will automatically rebuild and reload in the browser.

The static web site is output to the _site subdirectory. Do not commit it to the repository - we keep there only sources and neccessary configurations.

