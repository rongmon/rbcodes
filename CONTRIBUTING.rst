============
Contributing
============

Welcome to ``rbcodes`` contributor's guide.

This document focuses on getting any potential contributor familiarized
with the development processes, but `other kinds of contributions`_ are also
appreciated.

If you are new to using git_ or have never collaborated in a project previously,
please have a look at `contribution-guide.org`_. Other resources are also
listed in the excellent `guide created by FreeCodeCamp`_ [#contrib1]_.

Please notice, all users and contributors are expected to be **open,
considerate, reasonable, and respectful**. When in doubt, `Python Software
Foundation's Code of Conduct`_ is a good reference in terms of behavior
guidelines.


Issue Reports
=============

If you experience bugs or general issues with ``rbcodes``, please have a look
on the `issue tracker`_. If you don't see anything useful there, please feel
free to fire an issue report.

.. tip::
   Please don't forget to include the closed issues in your search.
   Sometimes a solution was already reported, and the problem is considered
   **solved**.

New issue reports should include information about your programming environment
(e.g., operating system, Python version) and steps to reproduce the problem.
Please try also to simplify the reproduction steps to a very minimal example
that still illustrates the problem you are facing. By removing other factors,
you help us to identify the root cause of the issue.


Documentation Improvements
==========================

You can help improve ``rbcodes`` docs by making them more readable and coherent, or
by adding missing information and correcting mistakes.

``rbcodes`` documentation is written in Markdown and kept in the ``docs/`` directory.
Documentation updates are contributed the same way as code contributions.

.. tip::
   Please notice that the `GitHub web interface`_ provides a quick way to
   propose changes in ``rbcodes``'s files. While this mechanism can
   be tricky for normal code contributions, it works perfectly fine for
   contributing to the docs, and can be quite handy.

   If you are interested in trying this method out, please navigate to
   the ``docs`` folder in the `source repository`_, find which file you
   would like to propose changes and click the little pencil icon at the
   top, to open `GitHub's code editor`_. Once you finish editing the file,
   please write a message in the form at the bottom of the page describing
   which changes you have made and what your motivations are,
   then submit your proposal.

When working on documentation changes in your local machine, you can
compile them using |tox|_::

    tox -e docs

and use Python's built-in web server for a preview in your web browser
(``http://localhost:8000``)::

    python3 -m http.server --directory 'docs/_build/html'


Code Contributions
==================

For more information about the project structure and main modules, please see the
`documentation <https://github.com/rongmon/rbcodes/blob/main/docs/main_readme.md>`_.

Submit an issue
---------------

Before you work on any non-trivial code contribution it's best to first create
a report in the `issue tracker`_ to start a discussion on the subject.
This often provides additional considerations and avoids unnecessary work.

Create an environment
---------------------

Before you start coding, we recommend creating an isolated `virtual
environment`_ to avoid any problems with your installed Python packages.
This can easily be done via Miniconda_::

    conda create -n rbcodes python=3.10
    conda activate rbcodes

Clone the repository
--------------------

#. Create a GitHub account if you do not already have one.
#. Fork the `rbcodes repository <https://github.com/rongmon/rbcodes>`_: click on the *Fork* button near the top of the
   page. This creates a copy of the code under your account on GitHub.
#. Clone your fork to your local disk::

    git clone https://github.com/<your-username>/rbcodes.git
    cd rbcodes

#. Install the package in development mode::

    pip install -U pip setuptools
    pip install -e .

   This allows you to import the package under development and test your changes.

Implement your changes
----------------------

#. Create a branch to hold your changes::

    git checkout -b my-feature

   Never work directly on the main branch!

#. Make your changes on this branch. Please add docstrings_ to new
   functions, modules and classes, especially if they are part of public APIs.

#. When you're done editing, commit your changes with a clear message::

    git add <MODIFIED FILES>
    git commit -m "Description of your changes"

   Writing a `descriptive commit message`_ is highly recommended.
   You can check the commit history with::

      git log --oneline -10

   to see the project's commit message style.

#. Add unit tests and documentation if your contribution adds a new feature
   (not just a bugfix).

#. Before submitting, check that your changes don't break unit tests::

    pytest

   or run the full test suite with::

    tox

   (install with ``pip install tox`` if needed).

Submit your contribution
------------------------

#. Push your local branch to your fork on GitHub::

    git push -u origin my-feature

#. Go to the GitHub page of your fork and click the "Compare & pull request" button
   to create a Pull Request.

#. In the PR description, explain:
   - What problem does this solve or what feature does it add?
   - How did you test this change?
   - Any related issues (reference them with ``#123``)

   You can also open the PR as a draft first and mark it as ready for review
   after any feedback from CI checks or code review.

#. We'll review your PR and provide feedback. Thank you for contributing!


Troubleshooting
---------------

The following tips can be used when facing problems to build or test the
package:

#. Make sure to fetch all the tags from the upstream repository_.
   The command ``git describe --abbrev=0 --tags`` should return the version you
   are expecting. If you are trying to run CI scripts in a fork repository,
   make sure to push all the tags.
   You can also try to remove all the egg files or the complete egg folder, i.e.,
   ``.eggs``, as well as the ``*.egg-info`` folders in the ``src`` folder or
   potentially in the root of your project.

#. Sometimes |tox|_ misses out when new dependencies are added, especially to
   ``setup.cfg`` and ``docs/requirements.txt``. If you find any problems with
   missing dependencies when running a command with |tox|_, try to recreate the
   ``tox`` environment using the ``-r`` flag. For example, instead of::

    tox -e docs

   Try running::

    tox -r -e docs

#. Make sure to have a reliable |tox|_ installation that uses the correct
   Python version (e.g., 3.7+). When in doubt you can run::

    tox --version
    # OR
    which tox

   If you have trouble and are seeing weird errors upon running |tox|_, you can
   also try to create a dedicated `virtual environment`_ with a |tox|_ binary
   freshly installed. For example::

    virtualenv .venv
    source .venv/bin/activate
    .venv/bin/pip install tox
    .venv/bin/tox -e all

#. `Pytest can drop you`_ in an interactive session in the case an error occurs.
   In order to do that you need to pass a ``--pdb`` option (for example by
   running ``tox -- -k <NAME OF THE FALLING TEST> --pdb``).
   You can also setup breakpoints manually instead of using the ``--pdb`` option.


Maintainer tasks
================

Releases
--------

For maintainers, releasing a new version involves:

#. Make sure all unit tests pass.
#. Tag the current commit on the main branch with a release tag, e.g., ``v2.0.0``.
#. Push the new tag to the repository::

    git tag v2.0.0
    git push origin v2.0.0

#. GitHub will automatically create a release from the tag. You can add release notes
   describing the changes in that release.

Automated releases to PyPI can be set up via GitHub Actions if desired in the future.



.. [#contrib1] Even though, these resources focus on open source projects and
   communities, the general ideas behind collaborating with other developers
   to collectively create software are general and can be applied to all sorts
   of environments, including private companies and proprietary code bases.


.. Links
.. _issue tracker: https://github.com/rongmon/rbcodes/issues
.. _source repository: https://github.com/rongmon/rbcodes


.. |virtualenv| replace:: ``virtualenv``
.. |pre-commit| replace:: ``pre-commit``
.. |tox| replace:: ``tox``


.. _black: https://pypi.org/project/black/
.. _CommonMark: https://commonmark.org/
.. _contribution-guide.org: https://www.contribution-guide.org/
.. _creating a PR: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request
.. _descriptive commit message: https://chris.beams.io/posts/git-commit
.. _docstrings: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
.. _first-contributions tutorial: https://github.com/firstcontributions/first-contributions
.. _flake8: https://flake8.pycqa.org/en/stable/
.. _git: https://git-scm.com
.. _GitHub's fork and pull request workflow: https://guides.github.com/activities/forking/
.. _guide created by FreeCodeCamp: https://github.com/FreeCodeCamp/how-to-contribute-to-open-source
.. _Miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _MyST: https://myst-parser.readthedocs.io/en/latest/syntax/syntax.html
.. _other kinds of contributions: https://opensource.guide/how-to-contribute
.. _pre-commit: https://pre-commit.com/
.. _PyPI: https://pypi.org/
.. _PyScaffold's contributor's guide: https://pyscaffold.org/en/stable/contributing.html
.. _Pytest can drop you: https://docs.pytest.org/en/stable/how-to/failures.html#using-python-library-pdb-with-pytest
.. _Python Software Foundation's Code of Conduct: https://www.python.org/psf/conduct/
.. _reStructuredText: https://www.sphinx-doc.org/en/master/usage/restructuredtext/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _tox: https://tox.wiki/en/stable/
.. _virtual environment: https://realpython.com/python-virtual-environments-a-primer/
.. _virtualenv: https://virtualenv.pypa.io/en/stable/

.. _GitHub web interface: https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files
.. _GitHub's code editor: https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files
