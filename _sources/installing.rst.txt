.. _installing:

Installation
============

PyPI
----

All third-party dependencies are installed automatically via pip, simply run:


.. code-block:: bash

    python -m pip install deltapd


Docker
------

The container must be able to access the input data, and write output data.
Therefore, you must specify a volume mapping (``-v host:container``), e.g.
the host path ``/tmp/deltapd`` is mapped to the container path ``/data``

Note: ``--plot`` will not work as it requires a display for PyQt.


.. code-block:: bash

    docker run -v /tmp/deltapd:/data deltapd:latest --r_tree /data/ar122_reference.tree --q_tree /data/ar122_ssu_query.tree --metadata /data/ar122_metadata_r95.tsv --msa_file /data/ar122_ssu_r95.fna --influence_thresh 2 --diff_thresh 0.10 --out_dir /data/deltapd_out --k 20 --cpus 10
