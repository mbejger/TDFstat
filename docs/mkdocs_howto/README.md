---
title: Documentation How To
nav_order: 1
---

# Documentation How To

## Introduction

This documentation is written in Markdown.
We use the [MkDocs](https://www.mkdocs.org/) static site generator (written in python) to create html pages.
The documentation sources (markdown) reside in the project's `docs` subdirectory (this is the mkdocs default).
Whenever you modify contents of this directory (and push changes to Github), Github runs our custom workflow to regenerate html pages and stores them in the gh-pages branch which is published at [https://polgraw.github.io/TDFstat/](https://polgraw.github.io/TDFstat/).

## Writing documentation

Some usefull documentation for authors:

- [Markdown](https://www.markdownguide.org/)
- [MkDocs documentation for users](https://www.mkdocs.org/user-guide/writing-your-docs/)
- [Math formulas](https://facelessuser.github.io/pymdown-extensions/extensions/arithmatex/) (`$` for inline and `$$` for block, `begin{equation}`, etc...)
- [Code highlighting](https://highlightjs.org/usage/) (e.g. ````python`)
- [Admonitions](https://python-markdown.github.io/extensions/admonition/)
  (`!!! note|danger|important|highlight|blink`)

Please [test any changes locally](#testing) before pushing them to Github!


## MkDocs configuration

The main configuration file for MkDocs, `mkdocs.yml`, resides in the project's root directory.

The nav section defines site menu.

In the theme section we enabled [Read the docs theme](https://www.mkdocs.org/user-guide/choosing-your-theme/#readthedocs).
RTD theme is modified using CSS file 'extra.css' and we change there the default width of content from 800 px to 100%. [Why the width is hardcoded in RTD?](https://github.com/readthedocs/sphinx_rtd_theme/issues/295).  
In addition we configure highlight.js package used by default in RTD. In particullar we load custom JS code from `javascripts/highlight.js`. This is the only way I have found to disable code highlighting by default in highlight.js, otherwise you get random languages detected on simple text. With this modification we have to specify the language explicitelly to enable highlighting (e.g. ````python`). Only [common languages are loaded by default](https://highlightjs.org/static/demo/), if you need [some other language](https://github.com/highlightjs/highlight.js/blob/main/SUPPORTED_LANGUAGES.md) add it to the list in `mkdocs.yml`.

We use MathJax along with Artihmatex extension (to preserve latex syntax with single $).
The whole configuration is build around new MathJax 3.x.
Following guidelines from [Artihmatex documentation](https://facelessuser.github.io/pymdown-extensions/extensions/arithmatex/#loading-mathjax) we have:

 - `generic` option is set to false  
 - we load tex-chtml combined components (for other options see [mathjax components docs](https://docs.mathjax.org/en/latest/web/components/combined.html))  
 - MathJax 3.x configuration is loaded via  `docs/javascripts/mathjax.js` file

!!! note
    There is no need to load MathJax explicitelly in mkdocs.yml since Arithmatex does this
!!! important
    `docs/javascripts/mathjax.js` must be loaded **before** MathJax (`- https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js`)
    

## Testing locally and pushing to Github {#testing}

You should install MkDocs in the conda virtual environment (with conda-forge enabled):
```
conda create -n mkdocs-test
conda activate mkdocs-test
conda install mkdocs pymdown-extensions pygments
```
Then go to the root directory of the repository, `<some dir>/TDFstat` and run
```
mkdocs serve -v
```
This should start a local webserver at address `http://127.0.0.1:8000/`. Check it in the web browser.
The pages will regenerate and reload autonatically if there are any changes.

After testing if the pages are working locally you can submit changes to master
```
git add <some files>]
git commit -a -m "[docs] something"
git push
```
and then push the generated site to gh-pages branch
```
mkdocs gh-deploy
```



