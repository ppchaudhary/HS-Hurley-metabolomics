# Input data

The feature tables are not included in this repository. They are available from the corresponding author on reasonable request, subject to institutional approval.

To run the pipeline, place these four files in this folder:

| File | Description |
|---|---|
| `Serum_data_HS.csv` | Serum feature table: column 1 = m/z, column 2 = putative name (blank if unannotated), then one column of intensities per sample (45 samples, `<n>-HS`) |
| `Metadata_serum.csv` | Two columns: Sample, Hurley Stage (1/2/3) |
| `Skin_data_HS.csv` | Skin feature table, same layout (15 samples, `<n>-MT`; intensities averaged per patient) |
| `Metadata_skin.csv` | Two columns: Sample, Hurley Stage |

Annotations are putative (accurate mass). Sample names in each feature table must match those in its metadata file.