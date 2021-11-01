.. selfies documentation master file, created by
   sphinx-quickstart on Sun Jun 14 23:40:28 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the SELFIES documentation!
=====================================

SELFIES (SELF-referencIng Embedded Strings) is a 100% robust molecular string
representation. A main objective is to use SELFIES as direct input into
machine learning models, in particular in generative models, for the
generation of outputs with guaranteed validity.

This library is intended to be light-weight and easy to use.

For explanation of the underlying principle (formal grammar) and experiments,
please see the `original paper`_.

For comments, bug reports or feature ideas, please use github issues
or send an email to mario.krenn@utoronto.ca and alan@aspuru.com.

.. _original paper: https://doi.org/10.1088/2632-2153/aba947

Installation
############

Install SELFIES in the command line using pip:

.. code-block::

   $ pip install selfies


.. toctree::
   :maxdepth: 2
   :caption: Contents

   derivation
   tutorial.ipynb
   selfies

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
