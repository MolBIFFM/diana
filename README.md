# pipeline

protein-protein interaction network assembly and analysis from mass spectrometry data

## setup
```
pip3 install -r pipeline/requirements.txt 
```
Developed using Python 3.9.7, Ubuntu 21.10 and Cytoscape 3.9.1.

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
A configuration file specifies a list of workflows executed sequentially. Configuration files are processed concurrently.

## configuration

The configuration specifies the assembly of a protein-protein interaction network from a set of input genes or proteins optionally associated with mass spectrometry data, using protein-protein interaction data from BioGRID, IntAct, MINT Reactome or STRING, optionally extended to proteins neighboring the input.

Enrichment of Gene Ontology terms or Reactome pathways by the protein-protein interaction network or its individual modules can be assessed, as well as the distribution of mass spectrometry measurements. 

Enrichment of proteins in the protein-protein interaction network by either Gene Ontology terms or Reactome pathways they are associated with can be assessed and exported with the network structure underlying these databases. Along with the networks, specific Cytoscape styles can be generated.

---

The specification of input genes or proteins.

```json
[
    {
      "genes": [
        {
          "file": null
        }
      ],
      "proteins": [
        {
          "file": null
        }
      ]
    }
]
```
The tabular input file.

```json
[
    {
      "genes": [
        {
          "accession column": null
        }
      ],
      "proteins": [
        {
          "accession column": null
        }
      ]
    }
]
```
The table column to extract UniProt gene or protein accessions from. These are mapped to primary UniProt accessions present in SwissProt. If they do not occur in SwissProt, they are discarded.

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
A regular expression used to extract all matching gene or protein accessions from an entry in the table, possibly stripping it of additional information. The default setting is `"^(.+?)$"`, corresponding to the entire entry.

```json
[
    {
      "genes": [
        {
          "sheet": 1
        }
      ],
      "proteins": [
        {
          "sheet": 1
        }
      ]
    }
]
```
The sheet of a spreadsheet to extract data from. The default setting is `1` corresponding to the first sheet of the file.

```json
[
    {
      "genes": [
        {
          "header": 1
        }
      ],
      "proteins": [
        {
          "header": 1
        }
      ]
    }
]
```
The line number of the header, allowing to skip lines. The default setting is `1`, corresponding to the first line of the sheet.

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

A list of input gene or protein accessions, alternative to extraction from tabular file.

```json
[
    {
      "genes": [
        {
          "taxonomy identifier": 9606
        }
      ]
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
The time of measurement to be associated with the measurements from the input file. The default setting is `0`.

```json
[
    {
      "proteins": [
        {
          "post-translational modification": "M"
        }
      ]
    }
]
```
An identifier for the type of post-translational modification associate with measurements from the corresponding file. The default setting is `"M"`.

```json
[
    {
      "proteins": [
        {
          "position column": null
        }
      ]
    }
]
```
The table column reporting modification sites of measurements. If available, measurements are sorted accordingly.

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
A regular expression used to extract all matching modification sites from an entry in the table, possibly stripping it of additional information. The default setting is `"^(.+?)$"`, corresponding to the entire entry.

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
A list of columns to extract replicate measurements from. The default setting is `[]`, corresponding to no data.

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
A threshold on the number of replicates required to consider a measurement. The default setting is `1`.

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
The maximum number of measurements to associate with each protein prioritized by largest absolute value. The default setting is `5`.

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
A function to combine individual replicates into a single measurement. The default setting is `"mean"`. Available settings are `"mean"` and `"median"`.

```json
[
    {
      "proteins": [
        {
          "logarithm": null
        }
      ]
    }
]
```
The base of the logarithm that measured measurements are expressed as. By default, ratios are assumed. Available settings are `null`, `2` and `10`.

```json
[
    {
      "networks": []
    }
]
```
A list of input protein-protein interaction networks as output by a workflow. If multiple networks are specified, they are combined.

---

The specification of sources of protein-protein interactions for the assembly of the protein-protein interaction network. Database-specific requirements can be defined, where each must be satisfied for an interaction to be incorporated.

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
An integer k specifying the extension of the network using proteins separated by up to k protein-protein interactions from the input proteins in the corresponding database. The default setting is 0, corresponding to no extension.

```json
[
    {
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
      }
    }
]
```
The NCBI taxonomy ID of the organism. The default and currently only fully supported setting is `9606`, corresponding to Homo sapiens.


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
        "BioGRID": {
          "version": null
        }
      }
    }
]
```
The version of the BioGRID database to use. The default setting is `null`, corresponding to the latest version.

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
A list of accepted PSI-MI terms for interaction types. The default setting is `[]`, corresponding to accepting any annotation.

```json
[
    {
      "protein-protein interactions": {
        "IntAct": {
          "score": 0.0
        },
        "MINT": {
          "score": 0.0
        }
      }
    }
]
```
A PSI-MI score threshold. The default setting is `0.0`.

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
          "neighborhood score": 0.0
        }
      }
    }
]
```
The STRING gene neighborhood score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "neighborhood transferred score": 0.0
        }
      }
    }
]
```
The STRING transferred gene neighborhood score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "fusion score": 0.0
        }
      }
    }
]
```
The STRING gene fusion score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "cooccurrence score": 0.0
        }
      }
    }
]
```
The STRING gene coooccurrence score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "coexpression score": 0.0
        }
      }
    }
]
```
The STRING gene coexpression score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "coexpression transferred score": 0.0
        }
      }
    }
]
```
The STRING transferred gene coexpression score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "experiments score": 0.0
        }
      }
    }
]
```
The STRING experiments score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "experiments transferred score": 0.0
        }
      }
    }
]
```
The STRING transferred experiments score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "database score": 0.0
        }
      }
    }
]
```
The STRING database score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "database transferred score": 0.0
        }
      }
    }
]
```
The STRING transferred database score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "textmining score": 0.0
        }
      }
    }
]
```
The STRING  textmining score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "textmining transferred score": 0.0
        }
      }
    }
]
```
The STRING transferred textmining score threshold. The default setting is `0.0`.

```json
[
    {
      "protein-protein interactions": {
        "STRING": {
          "combined score": 0.0
        }
      }
    }
]
```
The STRING combined score threshold. The default setting is `0.0`.

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
The version of STRING to use. The default setting is `11.5`.

---

The specification of Cytoscape styles. If not present, no Cytoscape styles are exported.

```json
[
    {
      "Cytoscape": {
        "bar chart": {
          "site combination": "absmax"
        },
        "node color": {
          "site combination": "absmax"
        }
      }
    }
]
```
The function used to derive protein-specific measurements from their individual sites. The default setting is `"absmax"`, corresponding to the largest absolute value. Available settings are `"mean"`, `"median"`, `"max"`, `"absmax"`, `"min"`, `"absmin"`, `"sum"`, `"abssum"` and `null`, such that sites are considered individually.

```json
[
    {
      "Cytoscape": {
        "bar chart": {
          "conversion": null
        },
        "node color": {
          "conversion": null
        }
      }
    }
]
```
The conversion of measurements that a range refers to. It defaults to the log2-fold measurement but may be set to `"standard score"` or `"quantile"` with respect to the distribution of a particular modification at a particular time of measurement across the protein-protein interaction network.

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
The range of the bar charts reporting measurements. The default setting is `[-1.0, 1.0]` if `"conversion"` is not set, `[-2.0, 2.0]` if `"conversion"` is set to `"standard score"` and `[0.025, 0.975]` if `"conversion"` is set to `"quantile"`.

```json
[
    {
      "Cytoscape": {
        "node color": {
          "measurement": [-1.0, 1.0]
        }
      },
    }
]
```
The range of combined measurements categorizing proteins by whether the range is exceeded or not. The default setting is `[-1.0, 1.0]` if `"conversion"` is not set, `[-2.0, 2.0]` if `"conversion"` is set to `"standard score"` and `[0.025, 0.975]` if `"conversion"` is set to `"quantile"`.

```json
[
    {
      "Cytoscape": {
        "edge transparency": null,
      }
    }
]
```
The function used to derive a combined edge confidence score from scores in IntAct, MINT and STRING. For lack of corresponding score, 1.0 is used for all interactions from BioGRID and Reactome. The combined score is reflected by edge transparency. By default any edge receives a score of 1.0. Available settings are `null`,  `"mean"`, `"median"`, `"max"`, `"min"`, `"sum"`, `"number"`. Additionally, `"BioGRID"`, `"IntAct"`, `"MINT"`, `"Reactome"` and `"STRING"` refer to the score in only the particular database.

---

The specification of Gene Ontology and Reactome enrichment analysis of the protein-protein interaction network.

```json
[
    {
      "Gene Ontology enrichment": {
        "test": "hypergeometric"
      },
      "Reactome enrichment": {
        "test": "hypergeometric"
      }
    }
]
```
The statistical test to assess enrichment. The default setting is `"hypergeometric"`.
Available settings are `"binomial"` and `"hypergeometric"`.

```json
[
    {
      "Gene Ontology enrichment": {
        "correction": "Benjamini-Hochberg"
      },
      "Reactome enrichment": {
        "correction": "Benjamini-Hochberg"
      }
    }
]
```
The procedure to correct p-values for multiple testing. The default setting is `"Benjamini-Hochberg"`.
Available settings are `"Benjamini-Hochberg"` and `"Bonferroni"`.

```json
[
    {
      "Gene Ontology enrichment": {
        "p": 1.0
      },
      "Reactome enrichment": {
        "p": 1.0
      }
    }
]
```
The corrected p-value threshold. The default setting is `1.0`.

```json
[
    {
      "Gene Ontology enrichment": {
        "taxonomy identifier": 9606
      },
       "Reactome enrichment": {
        "taxonomy identifier": 9606
      }
    }
]
```
The NCBI taxonomy ID of the organism. The default and currently only fully supported setting is `9606`, corresponding to Homo sapiens.

```json
[
    {
      "Gene Ontology enrichment": {
        "namespaces": [
          "cellular_component",
          "molecular_function",
          "biological_process"
        ]
      }
    }
]
```
The Gene Ontology namespaces to consider. The default setting is `["cellular_component", "molecular_function" "biological_process"]`.

---

The specification of a modular decomposition of the protein-protein interaction network using a network topological community detection algorithm.

```json
[
    {
      "module detection": {
        "module size": null
      }
    }
]
```
An upper bound on the number of nodes per module. Modules are iteratively subdivided until this threshold is met. The default setting is the number of proteins in the network, corresponding to a single iteration of the community detection algorithm.

```json
[
    {
      "module detection": {
        "module size combination": "mean"
      }
    }
]
```
The function to cobine sizes of modules into a value decisive to meeting the module size threshold. The default setting is `"mean"`. Available settings are `"mean"`, `"median"`, `"max"` and `"min"`.

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
The resolution parameter of modularity optimized by module detection. The default setting is `1.0`, corresponding to non-parameterized modularity. Larger values produce smaller modules, placing emphasis on the expected number of edges.

```json
[
    {
      "module detection": {
        "edge weight": null
      }
    }
]
```
The function used to derive a combined edge confidence score from scores in IntAct, MINT and STRING as well as 1.0 used for BioGRID and Reactome, respectively for a lack of corresponding scoring. The combined score is utilized as edge weight in module detection. By default any edge receives a score of 1.0, corresponding to an unweighted network. Available settings are `null`, `"mean"`, `"median"`, `"max"`, `"min"`, `"sum"`, `"number"`. Additionally, `"BioGRID"`, `"IntAct"`, `"MINT"`, `"Reactome"` and `"STRING"` refer to the score in only the particular database.

---

The specification of statistical tests on individual modules with respect to either Gene Ontology or Reactome enrichment or the distribution of measurements across the protein-protein interaction network. 

Gene Ontology and Reactome enrichment can be assessed with respect to the entire annotation or pathway map, respectively, or the protein-protein interaction network. 

To assess their distribution, measurements can be interpreted in a binary way measuring modules' enrichment of proteins which exhibit measurements exceeding a specified threshold. Alternatively, the distribution of measurements within separate modules can be compared with the remaining network.

The tests act as filter on the exported modules of the protein-protein interaction network. A module is only exported, if it is significant according to any of the specified tests.

```json
[
    {
      "module detection": {
        "Gene Ontology enrichment": {
          "annotation": {
            "test": "hypergeometric"
          },
          "network": {
            "test": "hypergeometric"
          }
        },
        "Reactome enrichment": {
          "map": {
            "test": "hypergeometric"
          },
          "network": {
            "test": "hypergeometric"
          }
        },
        "measurement enrichment": {
          "proteins": {
            "test": "hypergeometric"
          },
          "sites": {
            "test": "hypergeometric"
          }
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
        "measurement location": {
          "proteins": {
            "test": "Wilcoxon"
          },
          "sites": {
            "test": "Wilcoxon"
          }
        }
      }
    }
]
```
The statistical test to compare modification- and time-specific measurement distributions of each module in with the remaining network. The default and setting is `"Wilcoxon"`. Available settings are `"Welch"` and `"Wilcoxon"`.

```json
[
    {
      "module detection": {
        "Gene Ontology enrichment": {
          "annotation": {
            "correction": "Benjamini-Hochberg"
          },
          "network": {
            "correction": "Benjamini-Hochberg"
          }
        },
        "Reactome enrichment": {
          "map": {
            "correction": "Benjamini-Hochberg"
          },
          "network": {
            "correction": "Benjamini-Hochberg"
          }
        },
        "measurement enrichment": {
          "proteins": {
            "correction": "Benjamini-Hochberg"
          },
          "sites": {
            "correction": "Benjamini-Hochberg"
          }
        },
        "measurement location": {
          "proteins": {
            "correction": "Benjamini-Hochberg"
          },
          "sites": {
            "correction": "Benjamini-Hochberg"
          }
        }
      }
    }
]
```
The procedure to correct p-values for multiple testing. The default setting is `"Benjamini-Hochberg"`.
Available settings are `"Benjamini-Hochberg"` and `"Bonferroni"`.

```json
[
    {
      "module detection": {
        "Gene Ontology enrichment": {
          "annotation": {
            "p": 1.0
          },
          "network": {
            "p": 1.0
          }
        },
        "Reactome enrichment": {
          "map": {
            "p": 1.0
          },
          "network": {
            "p": 1.0
          }
        },
        "measurement enrichment": {
          "proteins": {
            "p": 1.0
          },
          "sites": {
            "p": 1.0
          }
        },
        "measurement location": {
          "proteins": {
            "p": 1.0
          },
          "sites": {
            "p": 1.0
          }
        }
      }
    }
]
```
The corrected p-value threshold. The default setting is `1.0`.

```json
[
    {
      "module detection": {
        "Gene Ontology enrichment": {
          "annotation": {
            "taxonomy identifier": 9606
          },
          "network": {
            "taxonomy identifier": 9606
          }
        },
        "Reactome enrichment": {
          "map": {
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
      }
    }
]
```
The Gene Ontology namespaces to consider. The default setting is `["cellular_component", "molecular_function" "biological_process"]`.

```json
[
    {
      "module detection": {
        "measurement enrichment": {
          "proteins": {
            "measurement": [-1.0, 1.0]
          },
          "sites": {
            "measurement": [-1.0, 1.0]
          }
        }
      }
    }
]
```
The range of measurements categorizing proteins by whether the range is exceeded or not. The default setting is `[-1.0, 1.0]` if `"conversion"` is not set, `[-2.0, 2.0]` if `"conversion"` is set to `"standard score"` and `[0.025, 0.975]` if `"conversion"` is set to `"quantile"`.

```json
[
    {
      "module detection": {
        "measurement enrichment": {
          "proteins": {
            "conversion": null
          },
          "sites": {
            "conversion": null
          }
        }
      }
    }
]
```
The conversion of measurements that a range refers to. It defaults to the log2-fold measurement but may be set to `"standard score"` or `"quantile"` with respect to the distribution of a particular modification at a particular time of measurement across the protein-protein interaction network.

```json
[
    {
      "module detection": {
        "measurement enrichment": {
          "proteins": {
            "site combination": "absmax"
          } 
        }
      }
    }
]
```
The function used to derive a protein-specific measurement from a its individual sites. The default setting is `"absmax"`, corresponding to the largest absolute value. Available settings are `"mean"`, `"median"`, `"max"`, `"absmax"`, `"min"`, `"absmin"`, `"sum"` and `"abssum"`.

---

The specification of Gene Ontology or Reactome network assembly from the protein-protein interaction network. 

A Gene Ontology network is composed of terms, a Reactome network of pathways represented protein-protein interaction network. They contain the enrichment of each respective entity by proteins from the protein-protein interaction network.

Optionally, the proteins considered may be restricted based on their associated measurements, either by a union or intersection of subsets of proteins exceeding the specified ranges.

```json
[
    {
      "Gene Ontology network": {
        "union": [
          {
            "time": null
          }
        ],
        "intersection": [
          {
            "time": null
          }
        ]
      },
      "Reactome network": {
        "union": [
          {
            "time": null
          }
        ],
        "intersection": [
          {
            "time": null
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
      "Gene Ontology network": {
        "union": [
          {
            "post-translational modification": null
          }
        ],
        "intersection": [
          {
            "post-translational modification": null
          }
        ]
      },
      "Reactome network": {
        "union": [
          {
            "post-translational modification": null
          }
        ],
        "intersection": [
          {
            "post-translational modification": null
          }
        ]
      }
    }
]
```
The modification considered to determine a subset of proteins.

```json
[
    {
      "Gene Ontology network": {
        "union": [
          {
            "site combination": "absmax"
          }
        ],
        "intersection": [
          {
            "site combination": "absmax",
          }
        ]
      },
      "Reactome network": {
        "union": [
          {
            "site combination": "absmax"
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
The function used to derive a protein-specific measurement from a its individual sites. The default setting is `"absmax"`, corresponding to the largest absolute value. Available settings are `"mean"`, `"median"`, `"max"`, `"absmax"`, `"min"`, `"absmin"`, `"sum"`, `"abssum"` and `null`, such that sites are considered individually.

```json
[
    {
      "Gene Ontology network": {
        "union": [
          {
            "conversion": null
          }
        ],
        "intersection": [
          {
            "conversion": null
          }
        ]
      },
      "Reactome network": {
        "union": [
          {
            "conversion": null
          }
        ],
        "intersection": [
          {
            "conversion": null
          }
        ]
      }
    }
]
```
The conversion of measurements that a range refers to. It defaults to the log2-fold measurement but may be set to `"standard score"` or `"quantile"` with respect to the distribution of a particular modification at a particular time of measurement across the protein-protein interaction network.

```json
[
    {
      "Gene Ontology network": {
        "union": [
          {
            "measurement": [-1.0, 1.0]
          }
        ],
        "intersection": [
          {
            "measurement": [-1.0, 1.0]
          }
        ]
      },
      "Reactome network": {
        "union": [
          {
            "measurement": [-1.0, 1.0]
          }
        ],
        "intersection": [
          {
            "measurement": [-1.0, 1.0]
          }
        ]
      }
    }
]
```
The range of combined measurements categorizing proteins by whether the range is exceeded or not. The default setting is `[-1.0, 1.0]` if `"conversion"` is not set, `[-2.0, 2.0]` if `"conversion"` is set to `"standard score"` and `[0.025, 0.975]` if `"conversion"` is set to `"quantile"`.

```json
[
    {
      "Gene Ontology network": {
        "test": "hypergeometric"
      },
      "Reactome network": {
        "test": "hypergeometric"
      }
    }
]
```
The statistical test to assess enrichment. The default setting is `"hypergeometric"`.
Available settings are `"binomial"` and `"hypergeometric"`.

```json
[
    {
      "Gene Ontology network": {
        "correction": "Benjamini-Hochberg"
      },
      "Reactome network": {
        "correction": "Benjamini-Hochberg"
      }
    }
]
```
The procedure to correct p-values for multiple testing. The default setting is `"Benjamini-Hochberg"`.
Available settings are `"Benjamini-Hochberg"` and `"Bonferroni"`.

```json
[
    {
      "Gene Ontology network": {
        "taxonomy identifier": 9606
      },
      "Reactome network": {
        "taxonomy identifier": 9606
      }
    }
]
```

The NCBI taxonomy ID of the organism. The default and currently only fully supported setting is `9606`, corresponding to Homo sapiens.

```json
[
    {
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
The Gene Ontology namespaces to consider. The default setting is `["cellular_component", "molecular_function" "biological_process"]`.

---
   
The configuration files provided refer to data sets supplemented with the following publications.

1. M. Hahn, A. Covarrubias-Pinto, L. Herhaus, S. Satpathy, K. Klann, K. B.
Boyle, C. Münch, K. Rajalingam, F. Randow, C. Choudhary, and I. Dikic, **SIK2 orchestrates actin-dependent host response upon *Salmonella* infection**, Proceedings of the National Academy of Sciences, vol. 118, no. 19, May 2021.

1. E. Fiskin, T. Bionda, I. Dikic, and C. Behrends, **Global Analysis of Host and Bacterial Ubiquitinome in Response to *Salmonella Typhimurium* Infection**, Molecular Cell, vol. 62, no. 6, pp. 967-981, Jun. 2016.

2. K. Klann, D. Bojkova, G. Tascher, S. Ciesek, C. Münch, and J. Cinatl, **Growth Factor Receptor Signaling Inhibition Prevents SARS-CoV-2 Replication**, Molecular Cell, vol. 80, no. 1, 164-174.e4, Oct. 2020.

3. C. Schmutz, E. Ahrné, C. A. Kasper, T. Tschon, I. Sorg, R. F. Dreier, A. Schmidt, and C. Arrieumerlou, **Systems-Level Overview of Host Protein Phosphorylation During *Shigella flexneri* Infection Revealed by Phosphoproteomics**, Molecular & Cellular Proteomics, vol. 12, no. 10, pp. 2952-2968, Oct. 2013.

---

The following resources can be accessed.

1. R. Oughtred, J. Rust, C. Chang, B. J. Breitkreutz, C. Stark, A. Willems, L. Boucher, G. Leung, N. Kolas, F. Zhang, S.Dolma, J. Coulombe-Huntington, A. Chatr-Aryamontri, K. Dolinski, M. Tyers, **The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions**, Protein Science, vol. 30, no. 1, pp. 187-200, Oct. 18
   
2. M. Ashburner, C. A. Ball, J. A. Blake, D. Botstein, H. Butler, J. M.Cherry, A. P. Davis, K. Dolinski, S. S. Dwight, J. T. Eppig, M. A. Harris, D. P. Hill, L. Issel-Tarver, A. Kasarskis, S. Lewis, J. C. Matese, J.  E. Richardson, M. Ringwald, G. M. Rubin, G. Sherlock, **Gene Ontology: tool for the unification of biology**, Nature Genetics,  vol. 25, no. 1, pp. 25-29, May 2000

3. The Gene Ontology Consortium, **The Gene Ontology resource: enriching a GOld mine**, Nucleic Acids Research, vol. 49, no. D1, pp. D325-D334, Jan. 2021
  
4. S. Orchard, M. Ammari, B. Aranda, L. Breuza, L. Briganti, F. Broackes-Carter, N. H. Campbell, G. Chavali, C. Chen, N del-Toro , M. Duesbury, M. Dumousseau, E. Galeota, U. Hinz, M. Iannuccelli, S. Jagannathan, R. Jimenez, J. Khadake, A. Lagreid, L. Licata, R. C. Lovering, B. Meldal, A. N. Melidoni, M. Milagros, D. Peluso, L. Perfetto, P. Porras, A. Raghunath, S. Ricard-Blum, B. Roechert, A. Stutz, M. Tognolli, K. van Roey, G. Cesareni, H. Hermjakob, **The MIntAct project--IntAct as a common curation platform for 11 molecular interaction databases**, Nucleic Acids Research, vol. 42, no. D1, pp. D358-D363, Jan. 2014
   
5. L. Licata, L. Briganti, D. Peluso, L. Perfetto, M. Iannuccelli, E. Galeota, F. Sacco, A. Palma, A. P. Nardozza, E. Santonico, L. Castagnoli, G. Cesareni, **MINT, the molecular interaction database: 2012 update**, Nucleic Acids Research, vol. 40, no. D1, pp. D857-D861, Jan. 2012
   
6. M. Gillespie, B. Jassal, R. Stephan, M. Milacic, K. Rothfels, A. Senff-Ribeiro, J. Griss, C. Sevilla, L. Matthews, C. Gong, C. Deng, T. Varusai, E. Ragueneau, Y. Haider, B. May, V. Shamovsky, J. Weiser, T. Brunson, N. Sanati, L. Beckman, X. Shao, A. Fabregat, K. Sidiropoulos, J. Murillo, G. Viteri, J. Cook, S. Shorser, G. Bader, E. Demir, C. Sander, R. Haw, G. Wu, L. Stein, H. Hermjakob, P. D’Eustachio, **The reactome pathway knowledgebase 2022**, Nucleic Acids Research, vol. 50, no. D1, pp.  D687-D692, Jan. 2022
   
7. D. Szklarczyk, A. L. Gable, D. Lyon, A. Junge, S. Wyder, J. Huerta-Cepas, M. Simonovic, N. T. Doncheva, J. H. Morris, P. Bork, L. J. Jensen, C. von Mering, **STRING v11: protein–protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets**, Nucleic Acids Research, vol. 47, no. D1, pp. D607-D613, Jan. 2019
    
8. The UniProt Consortium, **UniProt: the universal protein knowledgebase in 2021**, Nucleic Acids Research, vol. 49, no. D1, pp. D480-D489, Jan 2021

---

The following external libraries are utilized.

1. A. A. Hagberg, D. A. Schult and P. J. Swart, **Exploring network structure, dynamics, and function using NetworkX**, Proceedings of the 7th Python in Science Conference, pp. 11-15, Aug. 2008
   
2. W. McKinney, **Data Structures for Statistical Computing in Python**, Proceedings of the 9th Python in Science Conference, pp. 56-61, Jun. 2010.

3. P. Virtanen, R. Gommers, T. E. Oliphant, M. Haberland, T. Reddy, D. Cournapeau, E. Burovski, P. Peterson, W. Weckesser, J. Bright, S. J. van der Walt, M. Brett, J. Wilson, K. J. Millman, N. Mayorov, A. R. J. Nelson, E. Jones, R. Kern, Eric L., C. J. Carey, İ. Polat, Y. Feng, E. W. Moore, J. VanderPlas, D. Laxalde, J. Perktold, R. Cimrman, I. Henriksen, E. A. Quintero, C. R. Harris, A. M. Archibald, A. H. Ribeiro, F. Pedregosa, P. van Mulbregt, SciPy 1.0 Contributors,  **SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python**, Nature Methods, vol. 17 no.3, pp.261-272, Feb. 2020