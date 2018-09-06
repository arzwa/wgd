Visualization module
********************

The visualization module allows both interactive visualization using bokeh, as
well as generating static image files. Below a screenshot of the interactive
interface is included:

.. image:: ath_cpa_viz.png

The interactive interface allows modification of key parameters, such as the
histogram bin-width and KDE bandwidth. You are strongly encouraged to observe
the effects of modifications in these parameters, as they may reveal
visualization artifacts. As one can see from the screenshot, it allows
overlaying multiple distributions, overlaying histograms and KDEs, and
dynamically hiding and showing of distributions (by clicking the entries in the
legend). Note that to run the interactive visualization, a bokeh server should
be running, which you can initiate with the following command::

    bokeh serve &

Note that ``bokeh`` should be installed automatically when installing ``wgd``.

Alternatively, the ``viz`` module also allows generating static images when the
``--interactive`` flag is not set.


.. automodule:: wgd.viz
    :members:
    :private-members:
    :special-members: __init__