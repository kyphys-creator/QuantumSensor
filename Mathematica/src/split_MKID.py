#!/usr/bin/env python3
"""Split QSENSv2_MKID.nb into per-Section notebooks, keep only TiN content."""

from split_common import split_notebook

split_notebook(
    detector="MKID",
    src_name="QSENSv2_MKID.nb",
    material_label="TiN",
    remove_labels={
        # 04_material: drop TES (Al) Subsection
        "04_material": ["TES"],
        # 06_response_defs: drop Al Subsubsection inside "defined functions"
        "06_response_defs": ["Al"],
        # 07_response_plot: drop Al Subsubsection inside "plot"
        "07_response_plot": ["Al"],
    },
)
