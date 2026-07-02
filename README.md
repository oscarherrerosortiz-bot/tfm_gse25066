# Predicción transcriptómica de respuesta patológica completa en GSE25066

Repositorio asociado a mi Trabajo Fin de Máster en Bioinformática y Análisis de Datos
Biomédicos. El proyecto estudia si la expresión transcriptómica tumoral permite predecir
la respuesta patológica completa (`pCR`) frente a enfermedad residual (`RD`) tras
quimioterapia neoadyuvante en cáncer de mama.

Los datos proceden del estudio público
[GSE25066](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25066) medido con la
plataforma Affymetrix Human Genome U133A Array
[GPL96](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL96).

## Objetivo

Comparar un modelo basal de covariables, un modelo transcriptómico Elastic Net a nivel de
sonda y un modelo combinado. La evaluación incluye validación cruzada anidada,
leave-one-cohort-out (LOCO) y análisis de estabilidad de las sondas seleccionadas.

## Estructura del repositorio

Lo que se incluye en el repositorio:

```text
tfm_gse25066/
├── README.md          este documento
└── scripts/           flujo numerado, de la comprobación inicial a la estabilidad de sondas
    ├── 00_initial_dataset_check.py
    ├── 01_prepare_analysis_data.py
    ├── 02_eda_qc.py
    ├── 03a_generate_predictive_splits.R
    ├── 03b_run_predictive_evaluation.R
    └── 03c_probe_stability_analysis.R
```

Los datos, los informes, las tablas de resultados y las figuras no se incluyen: los datos
deben añadirse manualmente (ver más abajo) y las salidas se generan al ejecutar los scripts.

## Datos

Los datos no se redistribuyen y deben obtenerse desde NCBI GEO. Antes de ejecutar, crea las
carpetas de datos dentro de la raíz del proyecto, de modo que la estructura quede así:

```text
tfm_gse25066/
├── README.md
├── scripts/
└── data/
    └── raw/
        └── GSE25066/
            ├── GSE25066_series_matrix.txt.gz
            └── GPL96.annot.gz
```

Todas las rutas que utilizan los scripts son relativas a la raíz del proyecto, por lo que la
ubicación de esa carpeta en el ordenador es indiferente mientras se respete esta estructura
interna y los scripts se ejecuten desde la raíz. El primer script puede descargar la series
matrix si no existe localmente, y el segundo puede obtener la anotación GPL96; la descarga
automática requiere conexión a internet.

## Requisitos

Análisis en Python (preparación y EDA/QC) y R (modelado). Dependencias principales:

- Python 3.11: numpy, pandas, scipy, scikit-learn, matplotlib, requests.
- R: glmnet, pROC, data.table, digest, yaml. La comprobación inicial utiliza además
  `limma` (Bioconductor).

## Ejecución

Desde la raíz del repositorio, en este orden:

```bash
python scripts/00_initial_dataset_check.py
python scripts/01_prepare_analysis_data.py
python scripts/02_eda_qc.py
Rscript scripts/03a_generate_predictive_splits.R
Rscript scripts/03b_run_predictive_evaluation.R
Rscript scripts/03c_probe_stability_analysis.R
```

Cada script genera sus salidas (metadatos procesados, informes, tablas y figuras) en
carpetas locales que no se versionan.

## Limitaciones

Estudio retrospectivo y multicohorte con desbalance entre pCR y RD. La plataforma es un
microarray a nivel de sonda con ambigüedades de anotación. La evaluación LOCO estima
transportabilidad entre cohortes del mismo recurso, no una validación clínica externa
independiente. Los resultados son exploratorios, no causales y no demuestran utilidad
clínica directa.

## Referencia

Hatzis C, et al. A genomic predictor of response and survival following
taxane-anthracycline chemotherapy for invasive breast cancer. *JAMA*.
2011;305(18):1873-1881. [doi:10.1001/jama.2011.593](https://doi.org/10.1001/jama.2011.593).
