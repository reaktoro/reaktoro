<img src='resources/icons/reaktoro.svg' width='100%'>

[![Linux Build](https://github.com/reaktoro/reaktoro/workflows/linux/badge.svg?branch=master)](https://github.com/reaktoro/reaktoro/actions?query=workflow%3Alinux)
[![OSX Build](https://github.com/reaktoro/reaktoro/workflows/osx/badge.svg?branch=master)](https://github.com/reaktoro/reaktoro/actions?query=workflow%3Aosx)
[![Windows Build](https://github.com/reaktoro/reaktoro/workflows/windows/badge.svg?branch=master)](https://github.com/reaktoro/reaktoro/actions?query=workflow%3Awindows)
[![image](https://anaconda.org/conda-forge/reaktoro/badges/installer/conda.svg)](https://conda.anaconda.org/conda-forge)
<!-- [![image](https://badge.fury.io/py/reaktoro.svg)](https://badge.fury.io/py/reaktoro) -->

# Introduction

Reaktoro is a unified framework for modeling chemically reactive systems.

Below are some features and modeling capabilities of Reaktoro:

* support to several [thermochemical databases](TODO_ADD_LINK_TO_SUPPORTED_THERMODYNAMIC_DATABASES);
  * PHREEQC
  * SUPCRT
  * SUPCRTBL
  * ThermoFun
* support to chemical equilibrium calculations with general constraints;
* support to chemical kinetic calculations;
* efficient numerical algorithms implemented using modern programming techniques;
* the chemical systems can contain any number of phases;
* no limitations on the number of chemical species in each phase;
* use of automatic differentiation for computation of derivatives with respect to virtually any variable or parameter.

> Note: Chemical kinetic calculations are not yet supported in v2.0rc. Use
> Reaktoro v1.0+ if you need them.

## Installation and Tutorials

For installation instructions, tutorials, and list of publications related to
this project, please access [reaktoro.org](http://www.reaktoro.org). This web
site describes how to download and install Reaktoro, and demonstrate some basic
usage.

## FAQ

### How do I ask a question about Reaktoro?

If you have questions about using or installing Reaktoro, please go to
[Reaktoro's GitHub
Issues](https://github.com/reaktoro/Reaktoro/issues/new) and let us
know. Please select the **question** label on the right side of the
issue pages. We'll do our best to answer your question as soon as
possible.

### How can I report a bug?

You got a bug and this is frustrating, we understand you. But don't
worry — we'll be happy to fix it for you (*provided it is indeed a
bug!*).

Before you report a bug, please check first if someone else has already
reported the same issue. If not, go to [Reaktoro's GitHub
Issues](https://github.com/reaktoro/Reaktoro/issues/new) and enter a
*descriptive title* and *write your issue with enough details*. Please
select the label **bug** on the right side of the page.

Please provide a [Minimum Reproducible
Example](https://stackoverflow.com/help/mcve%3E)? Please provide such an
example so that we can be more efficient in identifying the bug and
fixing it for you.

Have you heard about
[Markdown](https://guides.github.com/features/mastering-markdown/)?
Please use Markdown syntax when reporting your issues.

### How can I contribute to Reaktoro?

First, thanks for your interest in contributing to Reaktoro! You can do
so in many ways, from reporting bugs and writing tutorials to helping us
with code development. You might also consider **financially supporting
Reaktoro's development** by helping us extending the development team if
you plan to make Reaktoro an essential software component in your
company or academic group.

Read more on how to contribute to Reaktoro [here](CONTRIBUTING.rst).

## Contributors

You can see the list of awesome people who has contributed code to
Reaktoro in the [contributors
page](https://github.com/reaktoro/Reaktoro/graphs/contributors).

We would love to have you as a contributor too, see
[CONTRIBUTING](CONTRIBUTING.rst) for more information.

<!-- ## Developing Quick-Start  TODO: Update this section in README

In order to start developing, you'll need to build Reaktoro from
sources. There are two ways: install the dependencies manually, as
described [here](http://www.reaktoro.org/installation.html), or using
Conda. [Conda](https://conda.io/docs/) is a tool for managing packages,
dependencies and environments for multiple languages, including Python
and C++, and supporting multiple platforms: Windows, Linux and macOS. In
order to start developing Reaktoro using Conda, these are the steps:

1. Install Miniconda, pick the 64-bit installer that uses the latest Python
   version from: [conda.io/miniconda.html](https://conda.io/miniconda.html).
2. Add `conda-forge` as a channel:
   `conda config --append channels conda-forge`
3. Install `conda-devenv`: `conda install -n base conda-devenv`
4. Create an environment for Reaktoro, from the repository root
   directory: `conda devenv`
5. Activate the environment: `source activate reaktoro` from
   Linux/macOS or `activate reaktoro` from Windows
6. Create a `build` directory and call `cmake` from it (for now check
   the <span class="title-ref">.travis.yml</span> file for an example
   on CMake parameters), OR, on Windows, call the `inv msvc` task to
   generate a project under `build\msvc` directory, open it in the IDE
   and build the `INSTALL` project. (`inv` is short for `invoke`, from
   the [Invoke](https://www.pyinvoke.org/) tool.) -->

## Contributing

Reaktoro is an open-source project and **we need your help** to make it
even more successful.

You can contribute to the project in many ways:

* [reporting bugs](#reporting-bugs),
* [proposing new features](#proposing-new-features),
* [performing benchmark calculations](#performing-benchmark-calculations),
* [contributing with documentation and examples](#contributing-with-documentation-and-examples),
* [contributing with development](#contributing-with-development).

You can see the list of awesome people who have contributed to Reaktoro in the
[contributors page](https://github.com/reaktoro/reaktoro/graphs/contributors).

### Reporting bugs

You got a bug and this is frustrating, we understand you. But don't
worry — we'll be happy to fix it for you (*provided it is indeed a
bug!*).

To report a bug, please go to [Reaktoro's GitHub
Issues](https://github.com/reaktoro/reaktoro/issues/new) and enter a
*descriptive title* and *write your issue with enough details*.

Have you heard about a [Minimum Reproducible
Example](https://stackoverflow.com/help/minimal-reproducible-example)? Please
provide such an example so that we can be more efficient in identifying the bug
and fixing it for you.

Have you heard about
[Markdown](https://guides.github.com/features/mastering-markdown/)?
Please use Markdown syntax when writing your issues.

### Proposing new features

You have a wish list for the Reaktoro development team (e.g.,
implementing certain modeling capabilities, supporting some specific
thermodynamic databases) and you want to propose that to us. You can do
so by going to [Reaktoro's GitHub
Issues](https://github.com/reaktoro/reaktoro/issues/new) and describing
there your desired new feature.

We'll do our best to get your proposed new features implemented.
Understand, however, that we have limited resources and a tight schedule
for ongoing projects, so that your requested additions could take some
time to materialize depending on their complexity.

If you foresee Reaktoro could become an essential software component for
the scientific and engineering investigations of your company or
academic group, please consider **financially supporting its
development** — we are always open to new additions to the team!

### Performing benchmark calculations

You have used Reaktoro to perform some calculations for which you have
experimental data (or you have performed similar calculations using
other codes). **We would be very excited if you could share with us your
results!** And if the results you obtained with Reaktoro are not great,
please **make sure we know about this** and we'll be very happy to help
you making your calculations more accurate.

Get in touch with us by going to [Reaktoro's GitHub
Issues](https://github.com/reaktoro/reaktoro/issues/new) and telling us
what you want to share.

### Contributing with documentation and examples

You have used Reaktoro for a while and you want now to contribute with
some examples on how to use it for solving some specific problems. Thank
you for your interest. We appreciate your effort and willingness to
contribute.

Please go to [Reaktoro's GitHub
Issues](https://github.com/reaktoro/reaktoro/issues/new) and write a new issue,
detailing what you want to do.

We have opted to use
[Markdown](https://guides.github.com/features/mastering-markdown/) when writing
documentation.

### Contributing with development

Great, you have C++ and/or Python experience and you want to contribute
to Reaktoro with code development! We're very excited about your
decision in joining forces with us to develop new features, fix bugs,
improve performance, and other tasks you might want to explore.

Before you start working with your contribution, please let us know what
you want to do by going to [Reaktoro's GitHub
Issues](https://github.com/reaktoro/reaktoro/issues/new) and detailing
there your intended development contribution.

By discussing your idea first with us, this promotes a more efficient
workflow for both you and us when it comes the time to merge your
contribution.

We use the [Fork & Pull Request
workflow](https://gist.github.com/Chaser324/ce0505fbed06b947d962) to
ensure everything works smoothly. Below is summary of the necessary
steps you need to take:

1. Fork this repository
2. Clone the repository at your machine
3. Add your changes in a branch named after what's being done (`lower-case-with-hyphens`)
4. Make a pull request to `reaktoro/reaktoro`, targeting the `master` branch

That's all for now. Is there any thing missing here? If so, please [let us
know](https://github.com/reaktoro/reaktoro/issues/new).

## License

LGPL v2.1

Copyright (C) 2014-2021 Allan Leal

Reaktoro is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at
your option) any later version.

Reaktoro is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
License for more details.
