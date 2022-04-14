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

The configuration defines the assembly of a protein-protein interaction network from a set of input genes or proteins optionally associated with mass spectrometry data, using protein-protein interaction data from BioGRID, CORUM, IntAct, MINT, Reactome or STRING, optionally extended to proteins neighboring the input.

The distribution of mass spectrometry measurements and enrichment of Gene Ontology terms or Reactome pathways by the protein-protein interaction network or its modules can be assessed. The protein-protein interaction network as well as networks derived from the Gene Ontology and Reactome can be exported along comprehensive Cytoscape style specifications.

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
The table column to extract UniProt gene or protein accessions from. These are mapped to primary UniProt accessions or discarded if not present in SwissProt. Isoform identifiers are maintained on primary but not transferred from secondary accessions.

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
          "species": 9606
        }
      ]
    }
]
```
The NCBI taxonomy ID of the target species. The default and currently only completely supported setting is `9606`, corresponding to Homo sapiens.

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
An identifier for the type of post-translational modification associate with measurements from the corresponding file. The default setting is `"M"`. Up to two types of post-translational modification per time of measurement are supported in Cytoscape styles.

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
The table column reporting modification sites of measurements, according to which they are sorted. If an entry contains less positions than measurements, missing sites are substituted by 0. If an entry contains more positions than measurements, only as many as there are measurements are

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
          "replicate combination": "mean"
        }
      ]
    }
]
```
A function to combine individual replicates into a single measurement. The function is applied to ratios, not their log2. The default setting is `"mean"`. Available settings are `"mean"` and `"median"`.

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

The specification of sources of protein-protein interactions for the assembly of the protein-protein interaction network. The protein-protein interaction network is exported if any source is specified. Database-specific requirements can be defined, where each must be satisfied for an interaction to be incorporated.

```json
[
    {
      "protein-protein interactions": {
        "BioGRID": {
          "neighbors": 0
        },
        "CORUM": {
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
An integer k specifying the extension of the network using species-specific proteins which are separated by up to k protein-protein interactions from the input proteins in the corresponding database. The default setting is 0, corresponding to no extension.

```json
[
    {
      "protein-protein interactions": {
        "BioGRID": {
          "species": 9606
        },
        "CORUM": {
          "species": 9606
        },
        "IntAct": {
          "species": 9606
        }, 
        "MINT": {
          "species": 9606
        },
        "Reactome": {
          "species": 9606
        },
        "STRING": {
          "species": 9606
        }
      }
    }
]
```
The NCBI taxonomy ID of the target species. The default and currently only completely supported setting is `9606`, corresponding to Homo sapiens.


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
A list of accepted interaction throughput annotations. The default setting is `[]`, corresponding to any annotation.

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
A list of accepted experimental system annotations. The default setting is `[]`, corresponding to any annotation.

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
A list of accepted experimental system type annotations. The default setting is `[]`, corresponding to any annotation.

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
        "CORUM": {
          "purification methods": []
        }
      }
    }
]
```
A list of accepted PSI-MI identifiers or terms for protein complex purification methods. The default setting is `[]`, corresponding to any annotation.

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
A list of accepted PSI-MI identifiers or terms for interaction detection methods. The default setting is `[]`, corresponding to any annotation.

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
A list of accepted PSI-MI identifiers or terms for interaction types. The default setting is `[]`, corresponding to any annotation.

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
A list of accepted interaction context annotations. The default setting is `[]`, corresponding to any annotation.

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
A list of accepted interaction type annotations. The default setting is `[]`, corresponding to any annotation.

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
The function used to derive protein-specific measurements from their individual sites. The function is applied to ratios, not their log2. The default setting is `"absmax"`, corresponding to the largest absolute value. Available settings are `"mean"`, `"median"`, `"max"`, `"absmax"`, `"min"`, `"absmin"`, `"sum"`, `"abssum"` and `null`, such that sites are considered individually.

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
        "species": 9606
      },
       "Reactome enrichment": {
        "species": 9606
      }
    }
]
```
The NCBI taxonomy ID of the target species. The default and currently only completely supported setting is `9606`, corresponding to Homo sapiens.

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

The specification of the assembly of a Gene Ontology or Reactome network. A Gene Ontology network is composed of terms, a Reactome network of pathways. Both report the enrichment of each respective entity by proteins in the protein-protein interaction network along with their relations in these databases.

The proteins considered can be restricted based on their associated measurements, either by a union or intersection of subsets of proteins exceeding the specified ranges.

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
The function used to derive a protein-specific measurement from a its individual sites. The function is applied to ratios, not their log2. The default setting is `"absmax"`, corresponding to the largest absolute value. Available settings are `"mean"`, `"median"`, `"max"`, `"absmax"`, `"min"`, `"absmin"`, `"sum"`, `"abssum"` and `null`, such that sites are considered individually.

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
        "species": 9606
      },
      "Reactome network": {
        "species": 9606
      }
    }
]
```

The NCBI taxonomy ID of the target species. The default and currently only completely supported setting is `9606`, corresponding to Homo sapiens.

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

Gene Ontology and Reactome enrichment can be assessed with respect to the respective annotation or the protein-protein interaction network. 

To assess the distribution of measurements, these can be classified in a binary way to measure the modules' enrichment of proteins which exhibit measurements exceeding a specified threshold. Alternatively, the distribution of measurements within separate modules can be compared with the remaining network.

These tests act as filter on the exported modules of the protein-protein interaction network. A module is exported if it is significant according to any of the specified tests.

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
          "annotation": {
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
          "annotation": {
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
          "annotation": {
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
            "species": 9606
          },
          "network": {
            "species": 9606
          }
        },
        "Reactome enrichment": {
          "annotation": {
            "species": 9606
          },
          "network": {
            "species": 9606
          }
        }
      }
    }
]
```
The NCBI taxonomy ID of the target species. The default and currently only completely supported setting is `9606`, corresponding to Homo sapiens.

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
   
The configuration files provided refer to data sets supplemented with the following publications.

- Fiskin, E. et al. (2016) **Global Analysis of Host and Bacterial Ubiquitinome in Response to *Salmonella Typhimurium* Infection**, *Mol. Cell*, 62, 967-981.

- Hahn, M. et al. (2021) **SIK2 orchestrates actin-dependent host response upon *Salmonella* infection**, *Proc. Natl. Acad. Sci.*, 118.

- Klann K. et al. (2020) **Growth Factor Receptor Signaling Inhibition Prevents SARS-CoV-2 Replication**, *Mol. Cell*, 80, 164-174.

- Schmutz, C. et al. (2013) **Systems-Level Overview of Host Protein Phosphorylation During *Shigella flexneri* Infection Revealed by Phosphoproteomics**, *Mol. Cell. Proteom.*, 12, 2952-2968.

---

The following resources can be accessed.

- Ashburner, M. et al. (2000) **Gene Ontology: tool for the unification of biology**, *Nat. Genet.*, 25, 25-29.

- The Gene Ontology Consortium (2021) **The Gene Ontology resource: enriching a GOld mine**, *Nucleic Acids Res.*, 49, D325-D334.

- Gillespie, M. et al. (2022) **The reactome pathway knowledgebase 2022**, *Nucleic Acids Res.*, 50, D687-D692.
  
- Giurgiu, M et al. (2019) **CORUM: the comprehensive resource of mammalian protein complexes-2019**, *Nucleic Acids Res.*, 47, D559-D563
  
- Licata, L. et al. (2012) **MINT, the molecular interaction database: 2012 update**, *Nucleic Acids Res.*, 40, D857-D861.

- Orchard, S. et al. (2014) **The MIntAct project-IntAct as a common curation platform for 11 molecular interaction databases**, *Nucleic Acids Res.*, 42, D358-D363.
 
- Oughtred, R. et al. (2018) **The BioGRID database: A comprehensive biomedical resource of curated protein, genetic, and chemical interactions**, *Protein Sci.*, 30, 187-200.
   
- Szklarczyk, D. et al. (2019) **STRING v11: protein-protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets**, *Nucleic Acids Res.*, 47, D607-D613.
    
- The UniProt Consortium (2021) **UniProt: the universal protein knowledgebase in 2021**, *Nucleic Acids Res.*, 49, D480-D489.

---

The following applications are integrated.

- Shannon, P. et al. (2003) **Cytoscape: a software environment for integrated models of biomolecular interaction networks**, *Genome Res.*, 13, 2498-2504.

---

The following external libraries are utilized.

- Hagberg, A. A. et al. (2008) **Exploring network structure, dynamics, and function using NetworkX**, *Proceedings of the 7th Python in Science Conference*, 11-15.
   
- McKinney, W. (2010) **Data Structures for Statistical Computing in Python**, *Proceedings of the 9th Python in Science Conference*, 56-61.

- Virtanen, P. et al. (2020)  **SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python**, *Nat. Methods*, 17, 261-272.

---

Development was inspired by previous work combining the following applications.
 
- Maere, S. et al. (2005) **BiNGO: a Cytoscape plugin to assess overrepresentation of Gene Ontology categories in Biological Networks**, *Bioinformatics*, 21, 3448-3449.

- Morris, J. H. et al. (2011) **clusterMaker: a multi-algorithm clustering plugin for Cytoscape**, *BMC Bioinform.*, 12.
  
- Su, G. et al. (2010) **GLay: community structure analysis of biological networks**, *Bioinformatics*, 26, 3135-3137.
  

---

References for implemented algorithms are listed in the corresponding source code.
