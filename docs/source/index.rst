.. selfies documentation master file, created by
   sphinx-quickstart on Sun Jun 14 23:40:28 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SELFIES' documentation!
===================================

SELFIES (SELF-referencIng Embedded Strings) is a general-purpose,
sequence-based, robust representation of semantically constrained graphs.
It is based on a Chomsky type-2 grammar, augmented with two self-referencing
functions. A main objective is to use SELFIES as direct input into machine
learning models, in particular in generative models, for the generation of
outputs with high validity.

The code presented here is a concrete application of SELFIES in chemistry, for
the robust representation of molecules (see the original paper on
`arXiv <https://arxiv.org/abs/1905.13741>`_). This library is intended to
be light-weight and easy to use, and also, allow users the flexibility
to customize the SELFIES language to their needs.

For comments, bug reports or feature ideas, please send an email to
mario.krenn@utoronto.ca and alan@aspuru.com.

Installation
############

Install SELFIES in the command line using pip:

.. code-block::

   $ pip install selfies


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorial
   selfies_examples.ipynb
   selfies

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
