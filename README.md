# pipeline
protein-protein interaction network assembly and analysis

## command line interface
```
usage: pipeline.py [-h] -c CONFIGURATIONS [CONFIGURATIONS ...] [-l {CRITICAL,ERROR,WARNING,INFO,DEBUG}] [-p PROCESSES]

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIGURATIONS [CONFIGURATIONS ...], --configurations CONFIGURATIONS [CONFIGURATIONS ...]
                        configuration files
  -l {CRITICAL,ERROR,WARNING,INFO,DEBUG}, --level {CRITICAL,ERROR,WARNING,INFO,DEBUG}
                        logging level (default: INFO)
  -p PROCESSES, --processes PROCESSES
                        maximum number of concurrent processes (default: 4)
```
A configuration file specifies a list of workflows sequentially executed. Configuration files are processed concurrently.

## configuration

```json
[
    {
      "genes": [
        {
          "file": ""
        }
      ],
      "proteins": [
        {
          "file": ""
        }
      ]
    }
]
```
The file location of the input file.

```json
[
    {
      "genes": [
        {
          "accession column": ""
        }
      ],
      "proteins": [
        {
          "accession column": ""
        }
      ]
    }
]
```
The column reporting UniProt gene or protein accessions.

```json
[
    {
      "genes": [
        {
          "accession format": "^(.+?)$"
        }
      ],
      "proteins": [
        {
          "accession format": "^(.+?)$"
        }
      ]
    }
]
```
A regular expression used to obtain gene or protein accessions from a cell value in the table. The default setting is `"^(.+?)$"`, corresponding to the entire entry.

```json
[
    {
      "genes": [
        {
          "sheet": ""
        }
      ],
      "proteins": [
        {
          "sheet": ""
        }
      ]
    }
]
```
The sheet of a spreadsheet. The default setting is `0` corresponding to the first sheet.

```json
[
    {
      "genes": [
        {
          "header": 0
        }
      ],
      "proteins": [
        {
          "header": 0
        }
      ]
    }
]
```
The line number of the header to skip initial lines of the sheet. The default setting is `0`, corresponding to the first line.

```json
[
    {
      "genes": [
        {
          "accessions": []
        }
      ],
      "proteins": [
        {
          "accessions": []
        }
      ]
    }
]
```
A list of input gene or protein accessions.

```json
[
    {
      "genes": [
        {
          "taxonomy identifier": 9606
        }
      ],
      "protein-protein interactions": {
        "BioGRID": {
          "taxonomy identifier": 9606
        },
        "IntAct": {
          "taxonomy identifier": 9606
        }, 
        "MINT": {
          "taxonomy identifier": 9606
        },
        "Reactome": {
          "taxonomy identifier": 9606
        },
        "STRING": {
          "taxonomy identifier": 9606
        }
      },
      "Gene Ontology enrichment": {
        "taxonomy identifier": 9606
      },
      "module detection": {
        "Gene Ontology enrichment": {
          "annotation": {
            "taxonomy identifier": 9606
          },
          "network": {
            "taxonomy identifier": 9606
          }
        }
      }
    }
]
```
The NCBI taxonomy ID of the organism. The default and currently only fully supported setting is `9606`, corresponding to Homo sapiens.

```json
[
    {
      "proteins": [
        {
          "time": 0
        }
      ]
    }
]
```
The time of measurement to be associated with the changes from the corresponding file. The default setting is `0`.

```json
[
    {
      "proteins": [
        {
          "post-translational modification": ""
        }
      ]
    }
]
```
An identifier for the type of post-translational modification associate with changes from the corresponding file. The default setting is `""`.

```json
[
    {
      "proteins": [
        {
          "position column": ""
        }
      ]
    }
]
```
The table column corresponding to modification sites used, if available, to order the changes. The default setting is `""`.

```json
[
    {
      "proteins": [
        {
          "position format": "^(.+?)$"
        }
      ]
    }
]
```
A regular expression used to obtain modification sites from a corresponding entry in the table. The default setting is `"^(.+?)$"`, corresponding to the entire entry.

```json
[
    {
      "proteins": [
        {
          "replicate columns": []
        }
      ]
    }
]
```
A list of columns to parse modification changes from. The default setting is `[]`, corresponding to no data.

```json
[
    {
      "proteins": [
        {
          "sites": 5
        }
      ]
    }
]
```
The maximum number of modification sites to associate with each protein prioritized by largest absolute change. The default setting is `5`.

```json
[
    {
      "proteins": [
        {
          "replicates": 1
        }
      ]
    }
]
```
A threshold on the number of replicates to accept a measurement. The default setting is `1`.

```json
[
    {
      "proteins": [
        {
          "combine replicates": "mean"
        }
      ]
    }
]
```
A function to combine individual replicates into a single change. The default setting is `"mean"`. Available settings are `"mean"` and `"median"`.

```json
[
    {
      "proteins": [
        {
          "logarithm": 0
        }
      ]
    }
]
```
The base of the logarithm that changes are expressed as and converted from to log2 scale. By default, ratios are assumed.

```json
[
    {
      "networks": []
    }
]
```
A list of input protein-protein interaction networks in GraphML format.

```json
[
    {
      "protein-protein interactions": {
        "BioGRID": {
          "neighbors": 0
        },
        "IntAct": {
          "neighbors": 0
        }, 
        "MINT": {
          "neighbors": 0
        },
        "Reactome": {
          "neighbors": 0
        },
        "STRING": {
          "neighbors": 0
        }
      }
    }
]
```
An integer k specifying the extension of the network by proteins separated by up to k protein-protein interactions satisfying the requirements from the input proteins in the corresponding database. The default setting is 0, corresponding to no extension.


```json
[
    {
      "protein-protein interactions": {
        "BioGRID": {
          "interaction throughput": []
        }
      }
    }
]
```
A list of accepted interaction throughput annotations. The default setting is `[]`, corresponding to accepting any annotation.

```json
[
    {
      "protein-protein interactions": {
        "BioGRID": {
          "experimental system": []
        }
      }
    }
]
```
A list of accepted experimental system annotations. The default setting is `[]`, corresponding to accepting any annotation.

```json
[
    {
      "protein-protein interactions": {
        "BioGRID": {
          "experimental system type": []
        }
      }
    }
]
```
A list of accepted experimental system type annotations. The default setting is `[]`, corresponding to accepting any annotation.

```json
[
    {
      "protein-protein interactions": {
        "BioGRID": {
          "multi-validated physical": false
        }
      }
    }
]
```
If true, restrict query to multi-validated physical protein-protein interactions. The default setting is `false`.

```json
[
    {
      "protein-protein interactions": {
        "IntAct": {
          "interaction detection methods": []
        },
        "MINT": {
          "interaction detection methods": []
        }
      }
    }
]
```
A list of accepted PSI-MI terms for interaction detection methods. The default setting is `[]`, corresponding to accepting any annotation.

```json
[
    {
      "protein-protein interactions": {
        "IntAct": {
          "interaction types": []
        },
        "MINT": {
          "interaction types": []
        }
      }
    }
]
```
A list of accepted PSI-MI terms interaction types. The default setting is `[]`, corresponding to accepting any annotation.

```json
[
    {
      "protein-protein interactions": {
        "IntAct": {
          "MI score": []
        },
        "MINT": {
          "MI score": []
        }
      }
    }
]
```
The PSI-MI score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "Reactome": {
          "interaction context": []
        }
      }
    }
]
```
A list of accepted interaction context annotations. The default setting is `[]`, corresponding to accepting any annotation.

```json
[
    {
      "protein-protein interactions": {
        "Reactome": {
          "interaction type": []
        }
      }
    }
]
```
A list of accepted interaction type annotations. The default setting is `[]`, corresponding to accepting any annotation.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "neighborhood": 0.0
        }
      }
    }
]
```
The STRING gene neighborhood score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "neighborhood transferred": 0.0
        }
      }
    }
]
```
The STRING transferred gene neighborhood score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "fusion": 0.0
        }
      }
    }
]
```
The STRING gene fusion score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "cooccurrence": 0.0
        }
      }
    }
]
```
The STRING gene coooccurrence score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "coexpression": 0.0
        }
      }
    }
]
```
The STRING gene coexpression score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "coexpression transferred": 0.0
        }
      }
    }
]
```
The STRING transferred gene coexpression score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "experiments": 0.0
        }
      }
    }
]
```
The STRING experiments score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "experiments transferred": 0.0
        }
      }
    }
]
```
The STRING transferred experiments score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "database": 0.0
        }
      }
    }
]
```
The STRING database score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "database transferred": 0.0
        }
      }
    }
]
```
The STRING transferred database score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "textmining": 0.0
        }
      }
    }
]
```
The STRING  textmining score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "textmining transferred": 0.0
        }
      }
    }
]
```
The STRING transferred textmining score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "combined": 0.0
        }
      }
    }
]
```
The STRING combined score threshold. The default setting is `0.0`, corresponding to no threshold.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "physical": false
        }
      }
    }
]
```
If true, restrict query to physical protein-protein interactions. The default setting is `false`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "version": 11.5
        }
      }
    }
]
```
The version of STRING to query. The default setting is `11.5`.

```json
[
    {
      "Cytoscape": {
        "bar chart": {
          "site combination": "absmax",
        },
        "node color": {
          "site combination": "absmax",
        }
      },
      "module detection": {
        "change enrichment": {
          "site combination": "absmax",
        }
      },
      "Reactome network": {
        "union": [
          {
            "site combination": "absmax",
          }
        ],
        "intersection": [
          {
            "site combination": "absmax",
          }
        ]
      },
      "Gene Ontology network": {
        "union": [
          {
            "site combination": "absmax",
          }
        ],
        "intersection": [
          {
            "site combination": "absmax",
          }
        ]
      }
    }
]
```
The function used to derive a proteins' representative change from a its individual sites. The default setting is `"absmax"`, corresponding to the largest absolute change. Available settings are `"mean"`, `"median"`, `"max"`, `"absmax"`, `"min"`, `"absmin"`, `"sum"`, `"abssum"`.



```json
[
    {
      "Cytoscape": {
        "bar chart": {
          "conversion": ""
        },
        "node color": {
          "conversion": ""
        }
      },
      "module detection": {
        "change enrichment": {
          "conversion": ""
        }
      },
      "Reactome network": {
        "union": [
          {
            "conversion": ""
          }
        ],
        "intersection": [
          {
            "conversion": ""
          }
        ]
      },
      "Gene Ontology network": {
        "union": [
          {
            "conversion": ""
          }
        ],
        "intersection": [
          {
            "conversion": ""
          }
        ]
      }
    }
]
```
The conversion of changes that change threshold range refers to. It defaults to log2-fold change but may be set to `"standard score"` or `"quantile"` with respect to the distribution of a change of a particular modification at a particular time of measurement in the network.

```json
[
    {
      "Cytoscape": {
        "bar chart": {
          "range": [-1.0, 1.0],
        }
      }
    }
]
```
The range of the bar charts in Cytoscape. The default setting is `[-1.0, 1.0]` if `"conversion"` is not set, `[-2.0, 2.0]` if `"conversion"` is set to `"standard score"` and `[0.025, 0.975]` if `"conversion"` is set to `"quantile"`.

```json
[
    {
      "Cytoscape": {
        "node color": {
          "combined change": [-1.0, 1.0]
        }
      },
      "module detection": {
        "change enrichment": {
          "combined change": [-1.0, 1.0]
        }
      },
      "Reactome network": {
        "union": [
          {
            "combined change": [-1.0, 1.0]
          }
        ],
        "intersection": [
          {
            "combined change": [-1.0, 1.0]
          }
        ]
      },
      "Gene Ontology network": {
        "union": [
          {
            "combined change": [-1.0, 1.0]
          }
        ],
        "intersection": [
          {
            "combined change": [-1.0, 1.0]
          }
        ]
      }
    }
]
```
The range of combined changes distinguishing proteins by whether the range is exceeded or not. The default setting is `[-1.0, 1.0]` if `"conversion"` is not set, `[-2.0, 2.0]` if `"conversion"` is set to `"standard score"` and `[0.025, 0.975]` if `"conversion"` is set to `"quantile"`.

```json
[
    {
      "Cytoscape": {
        "edge transparency": "",
      }
    }
]
```
The function used to derive a combined edge confidence score from scores in IntAct, MINT and STRING as well as 1.0 used for BioGRID and Reactome, respectively for a lack of corresponding score. The combined score is reflected by edge transparency. By default any edge receives a score of 1.0. Available settings are `"mean"`, `"median"`, `"max"`, `"min"`, `"sum"`, `"number"`, `"BioGRID"`, `"IntAct"`, `"MINT"`, `"Reactome"`, `"STRING"`.

```json
[
    {
      "module detection": {
        "edge weight": ""
      }
    }
]
```
The function used to derive a combined edge confidence score from scores in IntAct, MINT and STRING as well as 1.0 used for BioGRID and Reactome, respectively for a lack of corresponding scoring. The combined score utilized as edge weight in module detection. By default any edge receives a score of 1.0, corresponding to an unweighted network. Available settings are `"mean"`, `"median"`, `"max"`, `"min"`, `"sum"`, `"number"`, `"BioGRID"`, `"IntAct"`, `"MINT"`, `"Reactome"`, `"STRING"`.

```json
[
    {
      "module detection": {
        "module size": 0
      }
    }
]
```
An upper bound on the number of nodes in any module. Modules are iteratively subdivided until this threshold is met. The default setting is the number of proteins in the network, corresponding to a single iteration.

```json
[
    {
      "module detection": {
        "module size combination": "mean"
      }
    }
]
```
The function of the sizes of modules decisive to meeting the module size threshold. The default setting is `"mean"`, the mean of all module sizes. Available settings are `"mean"`, `"median"`, `"max"` and `"min"`.

```json
[
    {
      "module detection": {
        "algorithm": "Louvain"
      }
    }
]
```
The community detection algorithm. The default setting is `"Louvain"`. Available settings are `"Clauset-Newman-Moore"` and `"Louvain"`.

```json
[
    {
      "module detection": {
        "resolution": 1.0
      }
    }
]
```
The resolution parameter of modularity optimized by module detection. The default setting is `1.0`, corresponding to non-parameterized modularity.

```json
[
    {
      "Gene Ontology enrichment": {
        "test": "hypergeometric"
      },
      "module detection": {
        "Gene Ontology enrichment": {
          "annotation": {
            "test": "hypergeometric"
          },
          "network": {
            "test": "hypergeometric"
          }
        },
        "change enrichment": {
          "test": "hypergeometric"
        }
      }
    }
]
```
The statistical test to assess enrichment. The default setting is `"hypergeometric"`.
Available settings are `"binomial"` and `"hypergeometric"`.

```json
[
    {
      "module detection": {
        "change tendency": {
          "test": "Wilcoxon"
        }
      }
    }
]
```
The statistical test to assess the changes of each module in relation to the entire network. The default and currently only available setting is `"Wilcoxon"`, the Wilcoxon rank sum test.

```json
[
    {
      "Gene Ontology enrichment": {
        "correction": "Benjamini-Hochberg"
      },
      "module detection": {
        "Gene Ontology enrichment": {
          "annotation": {
            "correction": "Benjamini-Hochberg"
          },
          "network": {
            "correction": "Benjamini-Hochberg"
          }
        },
        "change enrichment": {
          "correction": "Benjamini-Hochberg"
        },
        "change tendency": {
          "correction": "Benjamini-Hochberg"
        }
      }
    }
]
```
The procedure to adjust p-values for multiple testing. The default setting is `"Benjamini-Hochberg"`.
Available settings are `"Benjamini-Hochberg"` and `"Bonferroni"`.

```json
[
    {
      "Gene Ontology enrichment": {
        "p": 1.0
      },
      "module detection": {
        "Gene Ontology enrichment": {
          "annotation": {
            "p": 1.0
          },
          "network": {
            "p": 1.0
          }
        },
        "change enrichment": {
          "p": 1.0
        },
        "change tendency": {
          "p": 1.0
        }
      }
    }
]
```
The threshold for corrected p-values. The default setting is `1.0`.

```json
[
    {
      "Gene Ontology enrichment": {
        "namespaces": [
          "cellular_component",
          "molecular_function",
          "biological_process"
        ]
      },
      "module detection": {
        "Gene Ontology enrichment": {
          "annotation": {
            "namespaces": [
              "cellular_component",
              "molecular_function",
              "biological_process"
            ]
          },
          "network": {
            "namespaces": [
              "cellular_component",
              "molecular_function",
              "biological_process"
            ]
          }
        }
      },
      "Gene Ontology network": {
        "namespaces": [
            "cellular_component",
            "molecular_function",
            "biological_process"
        ]
      }
    }
]
```
The Gene Ontology namespaces to consider. The default setting is `["cellular_component", "molecular_function" "biological_process"]`, corresponding to the entire Gene Ontology.

```json
[
    {
      "Reactome network": {
        "union": [
          {
            "time": 0
          }
        ],
        "intersection": [
          {
            "time": 0
          }
        ]
      },
      "Gene Ontology network": {
        "union": [
          {
            "time": 0
          }
        ],
        "intersection": [
          {
            "time": 0
          }
        ]
      }
    }
]
```
The time of measurement considered to determine a subset of proteins.

```json
[
    {
      "Reactome network": {
        "union": [
          {
            "post-translational modification": ""
          }
        ],
        "intersection": [
          {
            "post-translational modification": ""
          }
        ]
      },
      "Gene Ontology network": {
        "union": [
          {
            "post-translational modification": ""
          }
        ],
        "intersection": [
          {
            "post-translational modification": ""
          }
        ]
      }
    }
]
```
The modification considered to determine a subset of proteins.

---

The configuration files provided in this repository refer to data sets supplemented with the following publications.

1. M. Hahn, A. Covarrubias-Pinto, L. Herhaus, S. Satpathy, K. Klann, K. B.
Boyle, C. Münch, K. Rajalingam, F. Randow, C. Choudhary, and I. Dikic, **SIK2 orchestrates actin-dependent host response upon *Salmonella* infection**, Proceedings of the National Academy of Sciences, vol. 118, no. 19, e2024144118, May 2021.

2. E. Fiskin, T. Bionda, I. Dikic, and C. Behrends, **Global Analysis of Host and Bacterial Ubiquitinome in Response to *Salmonella Typhimurium* Infection**, Molecular Cell, vol. 62, no. 6, pp. 967-981, Jun. 2016.

3. K. Klann, D. Bojkova, G. Tascher, S. Ciesek, C. Münch, and J. Cinatl, **Growth Factor Receptor Signaling Inhibition Prevents SARS-CoV-2 Replication**, Molecular Cell, vol. 80, no. 1, 164-174.e4, Oct. 2020.

4. C. Schmutz, E. Ahrné, C. A. Kasper, T. Tschon, I. Sorg, R. F. Dreier, A. Schmidt, and C. Arrieumerlou, **Systems-Level Overview of Host Protein Phosphorylation During *Shigella flexneri* Infection Revealed by Phosphoproteomics**, Molecular & Cellular Proteomics, vol. 12, no. 10, pp. 2952-2968, Oct. 2013.
