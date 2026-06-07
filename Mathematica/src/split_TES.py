#!/usr/bin/env python3
"""Split QSENSv2_TES.nb into per-Section notebooks, keep only Al content."""

from split_common import split_notebook

split_notebook(
    detector="TES",
    src_name="QSENSv2_TES.nb",
    material_label="Al",
    remove_labels={
        # 04_material: drop MKID (TiN) Subsection
        "04_material": ["MKID"],
        # 06_response: drop TiN Subsubsections (defined-functions + plot)
        "06_response": ["TiN"],
    },
)
