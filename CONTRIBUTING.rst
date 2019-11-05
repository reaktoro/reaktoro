Contributing
============

Reaktoro is an open-source project and **we need your help** to make it even more successful.

You can contribute to the project in many ways:

- by `reporting bugs <#reporting-bugs>`__,
- by `proposing new features <#proposing-new-features>`__,
- by `performing benchmark calculations <#performing-benchmark-calculations>`__,
- by `contributing with documentation and examples <#contributing-with-documentation-and-examples>`__,
- by `contributing with development <#contributing-with-development>`__.

You can see the list of awesome people who has contributed to Reaktoro in the `contributors page <https://github.com/reaktoro/Reaktoro/graphs/contributors>`__.

Reporting bugs
--------------

You got a bug and this is frustrating, we understand you. But don't worry — we'll be happy to fix it for you (*provided it is indeed a bug!*).

To report a bug, please go to `Reaktoro's GitHub Issues`_ and enter a *descriptive title* and *write your issue with enough details*. On the right side of the page, click on *Labels* and choose **bug**.

Have you heard about a `Minimum Reproducible Example`_? Please provide such an example so that we can be more efficient in identifying the bug and fixing it for you.

Have you heard about `Markdown`_? Please use Markdown syntax when writing your issues.

Proposing new features
----------------------

You have a wish list for the Reaktoro development team (e.g., implementing certain modeling capabilities, supporting some specific thermodynamic databases) and you want to propose that to us. You can do so by going to `Reaktoro's GitHub Issues`_ and describing there your desired new feature. Please choose the **enhancement** label for this issue.

We'll do our best to get your proposed new features implemented. Understand, however, that we have limited resources and a tight schedule for ongoing projects, so that your requested additions could take some time to materialize depending on their complexity.

If you foresee Reaktoro could become an essential software component for the scientific and engineering investigations of your company or academic group, please consider **financially supporting its development** — we are always open to new additions to the team!


Performing benchmark calculations
---------------------------------

You have used Reaktoro to perform some calculations for which you have experimental data (or you have performed similar calculations using other codes). **We would be very excited if you could share with us your results!** And if the results you obtained with Reaktoro are not great, please **make sure we know about this** and we'll be very happy to help you making your calculations more accurate.

Get in touch with us by going to `Reaktoro's GitHub Issues`_ and telling us what you want to share.

Contributing with documentation and examples
--------------------------------------------

You have used Reaktoro for a while and you want now to contribute with some examples on how to use it for solving some specific problems. Thank you for your interest. We appreciate your effort and willingness to contribute.

Please go to `Reaktoro's GitHub Issues`_ and write a new issue, detailing what you want to do. Please select the label **enhancement**.

We have opted to use `reStructuredText`_ when writing documentation. There are some formatting requirements you need to follow, such as not inserting line breaks within paragraphs, but relying instead in your text editor to wrap lines automatically.

.. TODO: We should have a dedicated document describing these formatting and other requirements and then just point here to that document.

Contributing with development
-----------------------------

Great, you have C++ and/or Python experience and you want to contribute to Reaktoro with code development! We're very excited about your decision in joining forces with us to develop new features, fix bugs, improve performance, and other tasks you might want to explore.

Before you start working with your contribution, please let us know what you want to do by going to `Reaktoro's GitHub Issues`_ and detailing there your intended development contribution. Please select the label **enhancement**.

By discussing your idea first with us, this promotes a more efficient workflow for both you and us when it comes the time to merge your contribution.

We use the `Fork & Pull Request workflow`_ to ensure everything works smoothly. Below is summary of the necessary steps you need to take:

1. Fork this repository
2. Clone the repository at your machine
3. Add your changes in a branch named after what's being done (``lower-case-with-hyphens``)
4. Make a pull request to ``reaktoro/Reaktoro``, targeting the ``master`` branch

But what about how you can do your changes? In order to start developing, you'll need to build Reaktoro from sources. There
are two ways: install the dependencies manually, as described `here
<http://www.reaktoro.org/installation.html>`_, or using Conda. `Conda
<https://conda.io/docs/>`_ is a tool for managing packages, dependencies and
environments for multiple languages, including Python and C++, and supporting
multiple platforms: Windows, Linux and macOS. In order to start developing
Reaktoro using Conda, these are the steps:

#. Install Miniconda, pick the 64-bit installer that uses the latest Python version from: `conda.io/miniconda.html <https://conda.io/miniconda.html>`_.
#. Add ``conda-forge`` as a channel: ``conda config --append channels conda-forge``
#. Install ``conda-devenv``: ``conda install -n base conda-devenv``
#. Create an environment for Reaktoro, from the repository root directory: ``conda devenv``
#. Activate the environment: ``source activate reaktoro`` from Linux/macOS or ``activate reaktoro`` from Windows
#. Install pre-commit in order to activate git hooks, checkers and formatters: ``pre-commit install``
#. Create a ``build`` directory and call ``cmake`` from it (for now check the `.travis.yml` file for an example on CMake parameters), OR, on Windows, call the ``inv msvc`` task to generate a project under ``build\msvc`` directory, open it in the IDE and build the ``INSTALL`` project. (``inv`` is short for ``invoke``, from the `Invoke <https://www.pyinvoke.org/>`_ tool.)

We follow some standards for the sake of code style, uniformity and quality. Don't worry, everything is provided through
``pre-commit``. If you choose to develop Reaktoro using Conda, then the steps depicted above include ``pre-commit``.
If you installed the dependecies manually, you have to install ``pre-commit`` just after you cloned Reaktoro's
repository. Just do in your console: ``pre-commit install``, then checkers and formatters are enable to verify every
modification you introduced in a commit automatically. If something is not according to the standards, then the commit
will fail, but ``pre-commit`` will modify the necessary part in order to fit it in our standards, then just commit
again and everything will be fine.

That's all for now. Is there any thing missing here? If so, please `let us know`__.

.. _Reaktoro's GitHub Issues: https://github.com/reaktoro/Reaktoro/issues/new
.. _Minimum Reproducible Example: https://stackoverflow.com/help/mcve>
.. _Markdown: https://guides.github.com/features/mastering-markdown/
.. _reStructuredText: http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
.. _Fork & Pull Request workflow: https://gist.github.com/Chaser324/ce0505fbed06b947d962

__ `Reaktoro's GitHub Issues`_
