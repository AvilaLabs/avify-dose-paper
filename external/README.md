# External inputs

## Reference dosimetry archive
Goorley, Kiger and Zamenhof (2002), Medical Physics 29, 145–156.
https://doi.org/10.1118/1.1428758
Obtain the supplementary archive through the publisher's article page.
Expected filename: supplementary_material_1_1428758-sup-0001.zip
SHA-256: 511aa05c4253c77b2778a2bd7aa52334a1fd996cca3d08b7c1899f32db38fd62

The publisher archive itself is not redistributed. Comparison values in controls/w006_reference.json cite archive members. To verify against your own exact copy:

```bash
W006_ARCHIVE=/absolute/path/to/supplementary_material_1_1428758-sup-0001.zip python3 scripts/verify.py
```

The extraction script requires LibreOffice and regenerates controls/w006_reference.json; use a separate working copy if regenerating it. A different archive hash is not the frozen reference.

## Transport nuclear data
The original simulations used OpenMC 0.15.3 with ENDF/B-VII.1 HDF5 data, including neutron, photon and thermal scattering data. Upstream resources: https://openmc.org/data/ and https://docs.openmc.org/en/stable/usersguide/data.html . The unchanged harness expects ~/nuclear-data/endfb-vii.1-hdf5/cross_sections.xml. OPENMC_CROSS_SECTIONS alone does not override that Python configuration. Install or link the library at that location in a dedicated environment.

No nuclear-data library is distributed here, and an exact file-by-file manifest of the original transport library is not included. Consequently this package reproduces analysis of recorded tallies; it does not promise bitwise recreation of the original Monte Carlo history stream from an arbitrary current library download.
