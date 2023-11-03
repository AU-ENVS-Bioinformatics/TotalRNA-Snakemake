# Changelog

## [1.1.0](https://www.github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/compare/v1.0.0...v1.1.0) (2023-11-03)


### Features

* add diamond support and Silva138 ([#13](https://www.github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/issues/13)) ([4272233](https://www.github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/commit/42722336a3d24eb6ddf0f004bcd75c718e3a017b))
* mark some files as temporaru ([99bcb6f](https://www.github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/commit/99bcb6fa19c9779cbe745cf0460fca19f2279dbb))


### Bug Fixes

* broken import ([b2e03d7](https://www.github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/commit/b2e03d762c38b76c312e0c24416a29688d9a9c4b))
* don't remove legacy dir ([3e5121c](https://www.github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/commit/3e5121c488f6e205cca3090a33e3196dd21c4e2f))
* keep temporary files after metarib ([23305ec](https://www.github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/commit/23305ecc851d7630cf512f2ed9b173c712202198))
* typo path ([12b0e69](https://www.github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/commit/12b0e69322c4fea01b114bdae564761bcf00fcb3))
* unexpected exit ([631453a](https://www.github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/commit/631453ac0007fbdfe56302d93e072f1e4ccb489b))
* warning for manual corrrection ([#15](https://www.github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/issues/15)) ([3afdcf2](https://www.github.com/AU-ENVS-Bioinformatics/TotalRNA-Snakemake/commit/3afdcf284b6c84c48b2e642274c9043b77484347))

## 1.0.1 (2023-02-20)

### Fixed

- new Trinity assembly path

## 1.0.0 (2023-01-10)


### Features

- Trim reads using trim-galore.
- Filtering SSU and LSU reads using sormerna and SILVA.
- Reconstructing ribosomal genes using Metarib.
- Checking the quality of the ribosomal assembly using QUAST.
- Mapping RNA contigs to reads using CoMW.
- Classifying reads taxonomically using BLAST, SILVA and CREST.
- Assembling non-rRNA reads, filtering noncoding RNA, mapping mRNA reads to contigs and aligning contigs to SWORD using CoMW.
