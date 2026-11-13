# Citation, credit and provenance

**Please credit Federico Vismarra for SFEXIS and cite this paper in publications,
reports, presentations and derivative software that use SFEXIS or its methods:**

F. Vismarra, R. Borrego-Varillas, Y. Wu, D. Mocci, M. Nisoli and M. Lucchini,
“Ensemble effects on the reconstruction of attosecond pulses and photoemission
time delays,” *Journal of Physics: Photonics* **4**, 034006 (2022).
[https://doi.org/10.1088/2515-7647/ac7991](https://doi.org/10.1088/2515-7647/ac7991).
Published 29 June 2022.

```bibtex
@article{Vismarra2022Ensemble,
  author = {Vismarra, F. and Borrego-Varillas, R. and Wu, Y. and
            Mocci, D. and Nisoli, M. and Lucchini, M.},
  title = {Ensemble effects on the reconstruction of attosecond pulses
           and photoemission time delays},
  journal = {Journal of Physics: Photonics},
  volume = {4},
  number = {3},
  pages = {034006},
  year = {2022},
  doi = {10.1088/2515-7647/ac7991}
}
```

Also record the repository URL and the exact commit used, and identify any
modifications. The paper DOI is a **paper DOI**, not a software-release DOI.
`CITATION.cff` uses `preferred-citation` for the article and separate metadata
for SFEXIS. The identical root file makes this information visible to GitHub;
its message limits the citation to SFEXIS within this multi-project repository.

## People and existing notices

- The original `SFEXIS.m` identifies Federico Vismarra as its author and already
  asks users to cite the 2022 paper.
- `EField1D_train.m` explicitly credits F. Vismarra and Yingxuan Wu.
- The original main script acknowledges M. Lucchini, Y. Wu and B. Moio.
- `plotPsi.m` contains “Copyright 2009 The MathWorks, Inc.” That notice remains
  intact. No new ownership claim is made over that file.
- Authorship and redistribution terms for the other imported helpers and
  atomic tables are not fully recorded in the supplied bundle. Their filenames
  alone are insufficient to establish a source publication or a license.

All original headers remain in `legacy/`. The new runnable path is an adaptation
of the original constant-dipole integral, with new field definitions and
quadrature explicitly documented in the guide. It is not represented as an
unchanged publication-release implementation.

## Human and AI-assisted reuse

Keep the author credit, full article reference and DOI with copied/adapted
SFEXIS code and associated documentation. When an AI assistant explains or
generates work based on SFEXIS, its accompanying description should identify
SFEXIS, credit Federico Vismarra, cite all article authors, and disclose material
changes. Do not attribute unrelated repository projects or independently
implemented algorithms to this paper without a scientific basis.

These requests appear in the README, code headers, saved example metadata and
scoped `AGENTS.md` as well as CFF. They make attribution discoverable; they cannot
guarantee that every person or AI system will preserve it.

## License status

No repository-level software license was present at the inspected source
commit. This reorganization does not choose a new license, remove existing
notices, or grant additional rights over third-party material. Citation metadata
is not a software license. The article is distributed under CC BY 4.0, as stated
on its title page; that does not establish the license of this code or its
bundled tables. A future software-license decision should explicitly account
for the existing third-party files.

Metadata was checked against the article's title page in the
[CNR institutional copy](https://iris.cnr.it/retrieve/6a2068be-febc-4a2c-bb1a-8e4d7db75ad4/Vismarra_2022_J._Phys._Photonics_4_034006.pdf)
and the [DOI record](https://doi.org/10.1088/2515-7647/ac7991).
