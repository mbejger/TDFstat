---
layout: default
title: Docs configuration
excerpt:
nav_order: 1
---

# Docs Configuration

This documentation is written in Markdown. We use the standard Github Pages pipeline to generate viewable web pages. The Github pipeline, by default uses [Jekyll](https://jekyllrb.com/) static site generator written in Rubby. TDFstat project on Github is configured to use the `/docs` directory in the master branch as  source to this pipeline.

In file `_config.yml` we use an external [Read The Docs Theme for Jekyll and GitHub Pages](https://github.com/carlosperate/jekyll-theme-rtd).\\
We customized the default theme's file `_layouts/base.html` as follows:
- added support for mathjax (with $ as delimiters)
```
  <script>
    MathJax = {
	tex: {
	    inlineMath: [ ['$','$'], ["\\(","\\)"] ],
	    displayMath: [ ['$$','$$'], ["\\[","\\]"] ],
	    /*processEscapes: false,*/
      }
    };
  </script>
  <script type="text/javascript" id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"> </script>
```
- changed the default RTD theme width from 800 px to 100% ( [why is it hardcoded in RTD?](https://github.com/readthedocs/sphinx_rtd_theme/issues/295) )
```
  <style>
    .wy-nav-content { max-width: none; }
  </style>
```

There is also the `Gemfile` file. It is the same as the default one used by github internally to build the site with jekyll. We add it here explicitelly becuse it is needed to build the site locally.

## Usefull documentation
- [Jekyll](https://jekyllrb.com/docs/)
- [Github pages](https://support.github.com/features/pages)
- [Theme documentation](https://www.embeddedlog.com/jekyll-theme-rtd)
