
.. mdinclude:: ../README.md


Main Functionalities
----------------------------

.. autosummary::
   :toctree:
   :nosignatures:
   :template: main_functionalities.rst
   :caption: Main functionalities

   ukbppp_dl.pgwas.keep_significant_qtls_from_region
   ukbppp_dl.pgwas.list_tar_files_in_region_folder
   ukbppp_dl.common.download_from_synapse

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Full API
---------------------------------------

.. toctree::
   :caption: Full API
   :maxdepth: 4
   :hidden:

   ukbppp_dl

.. autosummary::
   :recursive:
   :nosignatures:

   .. We could have just written "ukbppp_dl" here, but we want to skip this extra layer, especially since there are so few modules to write anyway

   ukbppp_dl.common
   ukbppp_dl.pgwas


Index
-------------------

* :ref:`genindex`