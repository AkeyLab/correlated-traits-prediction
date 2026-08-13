"""One module per published figure.

Each module exposes a ``main()`` that regenerates its outputs under
``results/``. They are driven by :mod:`correlated_traits.cli`; import them
lazily, because several pull in matplotlib.
"""
