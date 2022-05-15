# DIANA

DIANA is a command line tool for **D**ata **I**ntegration **A**nd
**N**etwork-based **A**nalysis for post-translational modification mass
spectrometry data.

The tool assembles protein-protein interaction networks from genes or proteins
associated with mass spectrometric measurements and, optionally, extends to
proteins neighboring the input, querying protein-protein interaction data from
BioGRID, CORUM, IntAct, MINT, Reactome, and STRING.

The enrichment of CORUM protein complexes, Gene Ontology terms, and Reactome
pathways by the protein-protein interaction network as well as the distribution
of mass spectrometric measurements across it, its communities and subsets of
proteins derived from mass spectrometric measurements can be assessed.

In addition to a protein-protein interaction network and its individual
communities, networks of Gene Ontology terms or Reactome pathways, capturing
their respective hierarchical dependencies and enrichment by the protein-protein
interaction network, can each be exported along corresponding Cytoscape style
specifications.

## Setup

```
pip3 install -r diana/requirements.txt
```
External dependencies consist of NetworkX, pandas, and SciPy. DIANA is currently
developed using Python 3.10.4, Ubuntu 22.04 and Cytoscape 3.9.1.

## Command Line Interface

```
python3 diana/diana.py --help
```

## Configuration

Input genes or proteins can be read from tabular input files.

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
The tabular input file with UniProt gene or protein accessions and, optionally,
additional mass spectrometry data.

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
The table column to extract UniProt gene or protein accessions from. These are
mapped to primary UniProt accessions or discarded if not present in Swiss-Prot.
Isoform identifiers are maintained on primary, but not transferred from
secondary accessions.

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
A regular expression used to extract all matching gene or protein accessions
from an entry in the table, possibly removing additional information. The
default setting is `"^(.+?)$"`, corresponding to the entire entry.

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
The sheet of a spreadsheet to extract data from. The default setting is `1`
corresponding to the first sheet of the file.

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
The line number of the header, allowing to skip lines. The default setting is
`1`, corresponding to the first line of the sheet.

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

A list of input gene or protein accessions, alternative to extraction from
tabular file.

```json
[
    {
      "networks": [
        {
          "network": null
        }
      ]
    }
]
```
A protein-protein interaction network as output by a workflow. The union of
input protein-protein interaction networks with proteins and genes is used.

```json
[
    {
      "genes": [
        {
          "organism": 9606
        }
      ],
      "proteins": [
        {
          "organism": 9606
        }
      ],
      "networks": [
        {
          "organism": 9606
        }
      ]
    }
]
```
The NCBI taxonomy ID of the organism of interest. The default and currently only
completely supported setting is `9606`, corresponding to Homo sapiens.

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
The time of measurement to be associated with the measurements from the input
file. The default setting is `0`.

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
An identifier for the type of post-translational modification associate with
measurements from the corresponding file. The default setting is `"M"`.
Currently, up to two types of post-translational modification per time of
measurement are supported in Cytoscape styles.

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
The table column reporting modification sites of measurements, according to
which they are sorted. If an entry contains less positions than measurements,
missing sites are substituted by 0. If an entry contains more positions than
measurements, only as many as there are measurements are used in order.

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
A regular expression used to extract all matching modification sites from an
entry in the table, possibly removing additional information. The default
setting is `"^(.+?)$"`, corresponding to the entire entry.

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
A list of columns to extract replicate measurements from. The default setting is
`[]`, corresponding to no association with mass spectrometry data.

```json
[
    {
      "proteins": [
        {
          "replicate format": "^(.+?)$"
        }
      ]
    }
]
```
A regular expression used to extract all matching replicate measurements from an
entry in the table, possibly removing additional information. The default
setting is `"^(.+?)$"`, corresponding to the entire entry.

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
A threshold on the number of replicates required to consider a measurement. The
default setting is `1`.

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
The maximum number of measurements to associate with each protein prioritized by
largest absolute value. The default setting is `5`.

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
A function to combine individual replicates into a single measurement. The
function is applied to ratios, not their log2. The default setting is `"mean"`.
Available settings are `"mean"`, `"median"`, `"max"`, `"maxabs"`, `"min"`,
`"minabs"`, `"sum"`, and `"sumabs"`.

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
The base of the logarithm that measured measurements are expressed as. By
default, ratios are assumed, corresponding to `null`. Available settings are
`null`, `2` and `10`.

---

The specification of sources of protein-protein interactions for the assembly of
the protein-protein interaction network. The protein-protein interaction network
is exported if any source is specified. Database-specific requirements can be
specified, where each must be satisfied for an interaction to be incorporated. A
protein-protein interaction network is only exported if any database is queried.

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
An integer specifying the extension of the network using species-specific
proteins which are separated by up to `"neighbors"` protein-protein interactions
from the input proteins in the corresponding database. The default setting is 0,
corresponding to no extension.

```json
[
    {
      "protein-protein interactions": {
        "BioGRID": {
          "organism": 9606
        },
        "CORUM": {
          "organism": 9606
        },
        "IntAct": {
          "organism": 9606
        },
        "MINT": {
          "organism": 9606
        },
        "Reactome": {
          "organism": 9606
        },
        "STRING": {
          "organism": 9606
        }
      }
    }
]
```
The NCBI taxonomy ID for the organism of interest. The default and currently
only completely supported setting is `9606`, corresponding to Homo sapiens.


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
A list of accepted interaction throughput annotations. The default setting is
`[]`, corresponding to any annotation.

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
A list of accepted experimental system annotations. The default setting is `[]`,
corresponding to any annotation.

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
A list of accepted experimental system type annotations. The default setting is
`[]`, corresponding to any annotation.

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
If `true`, restrict query to multi-validated physical protein-protein
interactions. The default setting is `false`.

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
The version of the BioGRID database to query. The default setting is `null`,
corresponding to the latest version.

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
A list of accepted PSI-MI identifiers or terms for protein complex purification
methods. The default setting is `[]`, corresponding to any annotation.

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
A list of accepted PSI-MI identifiers or terms for interaction detection
methods. The default setting is `[]`, corresponding to any annotation.

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
A list of accepted PSI-MI identifiers or terms for interaction types. The
default setting is `[]`, corresponding to any annotation.

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
A list of accepted interaction context annotations. The default setting is `[]`,
corresponding to any annotation.

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
A list of accepted interaction type annotations. The default setting is `[]`,
corresponding to any annotation.

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
The STRING transferred gene neighborhood score threshold. The default setting is
`0.0`.

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
The STRING gene cooccurrence score threshold. The default setting is `0.0`.

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
The STRING transferred gene coexpression score threshold. The default setting is
`0.0`.

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
The STRING transferred experiments score threshold. The default setting is
`0.0`.

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
The STRING textmining score threshold. The default setting is `0.0`.

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
If `true`, restrict query to physical protein-protein interactions. The default
setting is `false`.

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

---

The specification of Cytoscape styles. If not present, no Cytoscape styles are
exported.

```json
[
    {
      "Cytoscape": {
        "bar chart": {
          "site combination": "maxabs"
        },
        "node color": {
          "site combination": "maxabs"
        }
      }
    }
]
```
The function used to derive protein-specific measurements from their individual
sites. The function is applied to ratios, not their log2. The default setting is
`"maxabs"`, corresponding to the largest absolute value. Available settings are
`"mean"`, `"median"`, `"max"`, `"maxabs"`, `"min"`, `"minabs"`, `"sum"`,
`"sumabs"` as well as `"increase"` and `"decrease"`, referring to the proportion
of a proteins' sites exhibiting either, and `null`, so that sites are considered
individually.

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
The conversion of measurements that a range refers to. It defaults to the
log2-fold measurement but may be set to `"standard score"` or `"quantile"`,
computed with respect to the distribution of a particular modification at a
particular time of measurement across the protein-protein interaction network.

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
The range of the bar charts reporting measurements. The default setting is
`[-1.0, 1.0]` if `"conversion"` is not set, `[-2.0, 2.0]` if `"conversion"` is
set to `"standard score"` and `[0.025, 0.975]` if `"conversion"` is set to
`"quantile"`.

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
The range of combined measurements categorizing proteins by whether the range is
exceeded or not. The default setting is `[-1.0, 1.0]` if `"conversion"` is not
set, `[-2.0, 2.0]` if `"conversion"` is set to `"standard score"` and
`[0.025, 0.975]` if `"conversion"` is set to `"quantile"`.

```json
[
    {
      "Cytoscape": {
        "edge transparency": null,
      }
    }
]
```
The function used to derive a combined edge confidence score from scores in
IntAct, MINT and STRING. For lack of corresponding score, 1.0 is used for all
interactions from BioGRID, CORUM and Reactome. The combined score is reflected
by edge transparency. By default any edge receives a score of 1.0. Available
settings are `null`, `"mean"`, `"median"`, `"max"`, `"min"`, `"sum"`,
`"number"`, the number of queried databases supporting the interaction. Further,
`"BioGRID"`, `"CORUM"`, `"IntAct"`, `"MINT"`, `"Reactome"`, and `"STRING"` refer
to the score in the particular database.

---

CORUM protein complex, Gene Ontology term and Reactome pathway enrichment of the
protein-protein interaction network can be assessed.

The proteins considered can be restricted, based on mass spectrometric
associated measurements, either by a union or intersection of specified subsets
of proteins from the protein-protein interaction network.

```json
[
    {
      "CORUM enrichment": {
        "test": "hypergeometric"
      },
      "Gene Ontology enrichment": {
        "test": "hypergeometric"
      },
      "Reactome enrichment": {
        "test": "hypergeometric"
      }
    }
]
```
The statistical test to assess enrichment. The default setting is
`"hypergeometric"`. Available settings are `"binomial"` and `"hypergeometric"`.

```json
[
    {
      "CORUM enrichment": {
        "correction": "Benjamini-Hochberg"
      },
      "Gene Ontology enrichment": {
        "correction": "Benjamini-Hochberg"
      },
      "Reactome enrichment": {
        "correction": "Benjamini-Hochberg"
      }
    }
]
```
The procedure to correct p-values for multiple testing. The default setting is
`"Benjamini-Hochberg"`. Available settings are `"Benjamini-Hochberg"` and
`"Bonferroni"`.

```json
[
    {
      "CORUM enrichment": {
        "p": 1.0
      },
      "Gene Ontology enrichment": {
        "p": 1.0
      },
      "Reactome enrichment": {
        "p": 1.0
      }
    }
]
```
The corrected p-value threshold to report and output a community. The default
setting is `1.0`.

```json
[
    {
      "CORUM enrichment": {
        "organism": 9606
      },
      "Gene Ontology enrichment": {
        "organism": 9606
      },
       "Reactome enrichment": {
        "organism": 9606
      }
    }
]
```
The NCBI taxonomy ID of the organism of interest. The default and currently only
completely supported setting is `9606`, corresponding to Homo sapiens.

```json
[
    {
      "CORUM enrichment": {
        "subsets": [
          {
            "time": null
          }
        ]
      },
      "Gene Ontology enrichment": {
        "subsets": [
          {
            "time": null
          }
        ]
      },
      "Reactome enrichment": {
        "subsets": [
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
      "CORUM enrichment": {
        "subsets": [
          {
            "post-translational modification": null
          }
        ]
      },
      "Gene Ontology enrichment": {
        "subsets": [
          {
            "post-translational modification": null
          }
        ]
      },
      "Reactome enrichment": {
        "subsets": [
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
      "CORUM enrichment": {
        "subsets": [
          {
            "site combination": "maxabs"
          }
        ]
      },
      "Gene Ontology enrichment": {
        "subsets": [
          {
            "site combination": "maxabs"
          }
        ]
      },
      "Reactome enrichment": {
        "subsets": [
          {
            "site combination": "maxabs"
          }
        ]
      }
    }
]
```
The function used to derive a protein-specific measurement from a its individual
sites. The function is applied to ratios, not their log2. The default setting is
`"maxabs"`, corresponding to the largest absolute value. Available settings are
`"mean"`, `"median"`, `"max"`, `"maxabs"`, `"min"`, `"minabs"`, `"sum"`,
`"sumabs"` as well as `"increase"` and `"decrease"`, referring to the proportion
of a proteins' sites exhibiting either, and `null`, so that sites are considered
individually.

```json
[
    {
      "CORUM enrichment": {
        "subsets": [
          {
            "conversion": null
          }
        ]
      },
      "Gene Ontology enrichment": {
        "subsets": [
          {
            "conversion": null
          }
        ]
      },
      "Reactome enrichment": {
        "subsets": [
          {
            "conversion": null
          }
        ]
      }
    }
]
```
The conversion of measurements that a range refers to. It defaults to the
log2-fold measurement but may be set to `"standard score"` or `"quantile"`,
computed with respect to the distribution of a particular modification at a
particular time of measurement across the protein-protein interaction network.

```json
[
    {
      "CORUM enrichment": {
        "subsets": [
          {
            "measurement": [-1.0, 1.0]
          }
        ]
      },
      "Gene Ontology enrichment": {
        "subsets": [
          {
            "measurement": [-1.0, 1.0]
          }
        ]
      },
      "Reactome enrichment": {
        "subsets": [
          {
            "measurement": [-1.0, 1.0]
          }
        ]
      }
    }
]
```
The range of combined measurements categorizing proteins by whether the range is
exceeded or not. The default setting is `[-1.0, 1.0]` if `"conversion"` is not
set, `[-2.0, 2.0]` if `"conversion"` is set to `"standard score"` and
`[0.025, 0.975]` if `"conversion"` is set to `"quantile"`.

```json
[
    {
      "CORUM enrichment": {
        "intersection": false
      },
      "Gene Ontology enrichment": {
        "intersection": false
      },
      "Reactome enrichment": {
        "intersection": false
      }
    }
]
```
If `true`, compute enrichment with respect to the intersection of specified
subsets of proteins from the protein-protein interaction network instead of
their union. The default setting is `false`.

```json
[
    {
      "CORUM enrichment": {
        "purification methods": []
      }
    }
]
```
A list of accepted PSI-MI identifiers or terms for protein complex purification
methods. The default setting is `[]`, corresponding to any annotation.

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
The Gene Ontology namespaces to consider. The default setting is
`["cellular_component", "molecular_function" "biological_process"]`.

---

Networks of Gene Ontology terms or Reactome pathways can be assembled. Both
incorporate the enrichment of each respective entity by the protein-protein
interaction network with respect to the annotation specific to an organism of
interest along the respective hierarchical relations of entities in these
databases.

The proteins considered can be restricted, based on mass spectrometric
associated measurements, either by a union or intersection of specified subsets
of proteins from the protein-protein interaction network.

```json
[
    {
      "Gene Ontology network": {
        "subsets": [
          {
            "time": null
          }
        ]
      },
      "Reactome network": {
        "subsets": [
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
        "subsets": [
          {
            "post-translational modification": null
          }
        ]
      },
      "Reactome network": {
        "subsets": [
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
        "subsets": [
          {
            "site combination": "maxabs"
          }
        ]
      },
      "Reactome network": {
        "subsets": [
          {
            "site combination": "maxabs"
          }
        ]
      }
    }
]
```
The function used to derive a protein-specific measurement from a its individual
sites. The function is applied to ratios, not their log2. The default setting is
`"maxabs"`, corresponding to the largest absolute value. Available settings are
`"mean"`, `"median"`, `"max"`, `"maxabs"`, `"min"`, `"minabs"`, `"sum"`,
`"sumabs"` as well as `"increase"` and `"decrease"`, referring to the proportion
of a proteins' sites exhibiting either, and `null`, so that sites are considered
individually.

```json
[
    {
      "Gene Ontology network": {
        "subsets": [
          {
            "conversion": null
          }
        ]
      },
      "Reactome network": {
        "subsets": [
          {
            "conversion": null
          }
        ]
      }
    }
]
```
The conversion of measurements that a range refers to. It defaults to the
log2-fold measurement but may be set to `"standard score"` or `"quantile"`,
computed with respect to the distribution of a particular modification at a
particular time of measurement across the protein-protein interaction network.

```json
[
    {
      "Gene Ontology network": {
        "subsets": [
          {
            "measurement": [-1.0, 1.0]
          }
        ]
      },
      "Reactome network": {
        "subsets": [
          {
            "measurement": [-1.0, 1.0]
          }
        ]
      }
    }
]
```
The range of combined measurements categorizing proteins by whether the range is
exceeded or not. The default setting is `[-1.0, 1.0]` if `"conversion"` is not
set, `[-2.0, 2.0]` if `"conversion"` is set to `"standard score"` and
`[0.025, 0.975]` if `"conversion"` is set to `"quantile"`.

```json
[
    {
      "Gene Ontology network": {
        "intersection": false
      },
      "Reactome network": {
        "intersection": false
      }
    }
]
```
If `true`, compute enrichment with respect to the intersection of specified
subsets of proteins from the protein-protein interaction network instead of
their union. The default setting is `false`.

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
The statistical test to assess enrichment. The default setting is
`"hypergeometric"`. Available settings are `"binomial"` and `"hypergeometric"`.

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
The procedure to correct p-values for multiple testing. The default setting is
`"Benjamini-Hochberg"`. Available settings are `"Benjamini-Hochberg"` and
`"Bonferroni"`.

```json
[
    {
      "Gene Ontology network": {
        "organism": 9606
      },
      "Reactome network": {
        "organism": 9606
      }
    }
]
```

The NCBI taxonomy ID of the organism of interest. The default and currently only
completely supported setting is `9606`, corresponding to Homo sapiens.

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
The Gene Ontology namespaces to consider. The default setting is
`["cellular_component", "molecular_function" "biological_process"]`.

---

Communities of the protein-protein interaction network can be detected by
parameterized modularity maximization and iterative subdivision.

```json
[
    {
      "community detection": {
        "algorithm": "Louvain"
      }
    }
]
```
The community detection algorithm. The default setting is `"Louvain"`. Available
settings are `"Clauset-Newman-Moore"` and `"Louvain"`.

```json
[
    {
      "community detection": {
        "resolution": 1.0
      }
    }
]
```
The resolution parameter of modularity which is maximized. The default setting
is `1.0`. Larger resolutions generate smaller communities, emphasizing the
expected number of intra-community edges.

```json
[
    {
      "community detection": {
        "edge weight": null
      }
    }
]
```
The function used to derive a combined edge score from confidence scores in
IntAct, MINT, and STRING as well as 1.0 used for BioGRID, CORUM and Reactome,
respectively, for a lack of a corresponding score. The combined score is
utilized as edge weight in community detection. By default any edge receives a
score of 1.0, corresponding to an unweighted network. Available settings are
`null`, `"mean"`, `"median"`, `"max"`, `"min"`, `"sum"`, `"number"`, the number
of queried databases supporting the interaction. Further, `"BioGRID"`,
`"CORUM"`, `"IntAct"`, `"MINT"`, `"Reactome"`, and `"STRING"` refer to the score
in that particular database.

```json
[
    {
      "community detection": {
        "community size": null
      }
    }
]
```
An additional upper bound on the number of proteins per community. Modules are
iteratively subdivided until this threshold is met. The default setting is the
number of proteins in the network, resulting in a single iteration of the
community detection algorithm.

```json
[
    {
      "community detection": {
        "community size combination": "mean"
      }
    }
]
```
The function to combine sizes of communities into a value compared to the
community size threshold. The default setting is `"mean"`. Available settings
are `"mean"`, `"median"`, `"max"`, and `"min"`.

---

CORUM protein complex, Gene Ontology term and Reactome pathway enrichment of
individual communities can be assessed. The proteins considered can be
restricted, based on mass spectrometric associated measurements, either by a
union or intersection of specified subsets of proteins from the protein-protein
interaction network.

To assess the distribution of mass spectrometry measurements, these are either
categorized to measure the communities' enrichment of proteins which exhibit
measurements exceeding a specified absolute or relative threshold.
Alternatively, the distribution of measurements within separate communities can
be compared with the remaining network with respect to either proteins or
modification sites.

The statistical tests act as filter on the communities of the protein-protein
interaction network in that a community is exported only if it appears
significant with respect to any of the specified tests.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "test": "hypergeometric"
        },
        "Gene Ontology enrichment": {
          "test": "hypergeometric"
        },
        "Reactome enrichment": {
          "test": "hypergeometric"
        },
        "measurement enrichment": {
          "test": "hypergeometric"
        }
      }
    }
]
```
The statistical test to assess enrichment. The default setting is
`"hypergeometric"`. Available settings are `"binomial"` and `"hypergeometric"`.

```json
[
    {
      "community detection": {
        "measurement location": {
          "test": "Wilcoxon"
        }
      }
    }
]
```
The statistical test to compare modification- and time-specific measurement
distributions of each community in with the remaining network. The default and
setting is `"Wilcoxon"`. Available settings are `"Welch"` and `"Wilcoxon"`.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "correction": "Benjamini-Hochberg"
        },
        "Gene Ontology enrichment": {
          "correction": "Benjamini-Hochberg"
        },
        "Reactome enrichment": {
          "correction": "Benjamini-Hochberg"
        },
        "measurement enrichment": {
          "correction": "Benjamini-Hochberg"
        },
        "measurement location": {
          "correction": "Benjamini-Hochberg"
        }
      }
    }
]
```
The procedure to correct p-values for multiple testing. The default setting is
`"Benjamini-Hochberg"`. Available settings are `"Benjamini-Hochberg"` and
`"Bonferroni"`.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "p": 1.0
        },
        "Gene Ontology enrichment": {
          "p": 1.0
        },
        "Reactome enrichment": {
          "p": 1.0
        },
        "measurement enrichment": {
          "p": 1.0
        },
        "measurement location": {
          "p": 1.0
        }
      }
    }
]
```
The corrected p-value threshold. The default setting is `1.0`.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "organism": 9606
        },
        "Gene Ontology enrichment": {
          "organism": 9606
        },
        "Reactome enrichment": {
          "organism": 9606
        }
      }
    }
]
```
The NCBI taxonomy ID of the organism of interest. The default and currently only
completely supported setting is `9606`, corresponding to Homo sapiens.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "annotation": false
        },
        "Gene Ontology enrichment": {
          "annotation": false
        },
        "Reactome enrichment": {
          "annotation": false
        }
      }
    }
]
```
If `true`, compute enrichment with respect to the entire annotation, specific to
the organism of interest, otherwise with respect to the protein-protein
interaction network. The default setting is `false`.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "subsets": [
            {
              "time": null
            }
          ]
        },
        "Gene Ontology enrichment": {
          "subsets": [
            {
              "time": null
            }
          ]
        },
        "Reactome enrichment": {
          "subsets": [
            {
              "time": null
            }
          ]
        }
      }
    }
]
```
The time of measurement considered to determine a subset of proteins.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "subsets": [
            {
              "post-translational modification": null
            }
          ]
        },
        "Gene Ontology enrichment": {
          "subsets": [
            {
              "post-translational modification": null
            }
          ]
        },
        "Reactome enrichment": {
          "subsets": [
            {
              "post-translational modification": null
            }
          ]
        }
      }
    }
]
```
The modification considered to determine a subset of proteins.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "subsets": [
            {
              "site combination": "maxabs"
            }
          ]
        },
        "Gene Ontology enrichment": {
          "subsets": [
            {
              "site combination": "maxabs"
            }
          ]
        },
        "Reactome enrichment": {
          "subsets": [
            {
              "site combination": "maxabs"
            }
          ]
        }
      }
    }
]
```
The function used to derive a protein-specific measurement from a its individual
sites. The function is applied to ratios, not their log2. The default setting is
`"maxabs"`, corresponding to the largest absolute value. Available settings are
`"mean"`, `"median"`, `"max"`, `"maxabs"`, `"min"`, `"minabs"`, `"sum"`,
`"sumabs"` as well as `"increase"` and `"decrease"`, referring to the proportion
of a proteins' sites exhibiting either, and `null`, so that sites are considered
individually.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "subsets": [
            {
              "conversion": null
            }
          ]
        },
        "Gene Ontology enrichment": {
          "subsets": [
            {
              "conversion": null
            }
          ]
        },
        "Reactome enrichment": {
          "subsets": [
            {
              "conversion": null
            }
          ]
        }
      }
    }
]
```
The conversion of measurements that a range refers to. It defaults to the
log2-fold measurement but may be set to `"standard score"` or `"quantile"`,
computed with respect to the distribution of a particular modification at a
particular time of measurement across the protein-protein interaction network.
Conversions refer to individual communities.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "subsets": [
            {
              "measurement": [-1.0, 1.0]
            }
          ]
        },
        "Gene Ontology enrichment": {
          "subsets": [
            {
              "measurement": [-1.0, 1.0]
            }
          ]
        },
        "Reactome enrichment": {
          "subsets": [
            {
              "measurement": [-1.0, 1.0]
            }
          ]
        }
      }
    }
]
```
The range of combined measurements categorizing proteins by whether the range is
exceeded or not. The default setting is `[-1.0, 1.0]` if `"conversion"` is not
set, `[-2.0, 2.0]` if `"conversion"` is set to `"standard score"` and
`[0.025, 0.975]` if `"conversion"` is set to `"quantile"`.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "intersection": false
        },
        "Gene Ontology enrichment": {
          "intersection": false
        },
        "Reactome enrichment": {
          "intersection": false
        }
      }
    }
]
```
If `true`, compute enrichment with respect to the intersection of specified
subsets of proteins from the protein-protein interaction network instead of
their union. The default setting is `false`.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "purification methods": []
        }
      }
    }
]
```
A list of accepted PSI-MI identifiers or terms for protein complex purification
methods. The default setting is `[]`, corresponding to any annotation.

```json
[
    {
      "community detection": {
        "Gene Ontology enrichment": {
          "namespaces": [
              "cellular_component",
              "molecular_function",
              "biological_process"
            ]
        }
      }
    }
]
```
The Gene Ontology namespaces to consider. The default setting is
`["cellular_component", "molecular_function" "biological_process"]`.

```json
[
    {
      "community detection": {
        "measurement enrichment": {
          "measurement": [-1.0, 1.0]
        }
      }
    }
]
```
The range of measurements categorizing proteins by whether the range is exceeded
or not. The default setting is `[-1.0, 1.0]` if `"conversion"` is not set,
`[-2.0, 2.0]` if `"conversion"` is set to `"standard score"` and
`[0.025, 0.975]` if `"conversion"` is set to `"quantile"`.

```json
[
    {
      "community detection": {
        "measurement enrichment": {
          "conversion": null
        }
      }
    }
]
```
The conversion of measurements that a range refers to. It defaults to the
log2-fold measurement but may be set to `"standard score"` or `"quantile"`,
computed with respect to the distribution of a particular modification at a
particular time of measurement across the protein-protein interaction network.

```json
[
    {
      "community detection": {
        "measurement enrichment": {
          "proteins": {
            "site combination": "maxabs"
          }
        }
      }
    }
]
```
The function used to derive a protein-specific measurement from a its sites. The
default setting is `"maxabs"`, corresponding to the largest absolute value.
Available settings are `"mean"`, `"median"`, `"max"`, `"maxabs"`, `"min"`,
`"minabs"`, `"sum"` and `"sumabs"` as well as `"increase"` and `"decrease"`,
referring to the proportion of a proteins' sites exhibiting either, and `null`,
so that sites are considered individually.

## References

The configuration files in this repository refer to data sets supplemented with
the following publications.

- Fiskin, E. et al. (2016) **Global Analysis of Host and Bacterial**
  **Ubiquitinome in Response to *Salmonella* Typhimurium Infection**, *Mol.*
  *Cell*, 62, 967-981.

- Hahn, M. et al. (2021) **SIK2 orchestrates actin-dependent host response**
  **upon *Salmonella* infection**, *Proc. Natl. Acad. Sci.*, 118.

- Klann K. et al. (2020) **Growth Factor Receptor Signaling Inhibition**
  **Prevents SARS-CoV-2 Replication**, *Mol. Cell*, 80, 164-174.

- Schmutz, C. et al. (2013) **Systems-Level Overview of Host Protein**
  **Phosphorylation During *Shigella flexneri* Infection Revealed by**
  **Phosphoproteomics**, *Mol. Cell. Proteom.*, 12, 2952-2968.

---

The following resources can be accessed.

- Ashburner, M. et al. (2000) **Gene Ontology: tool for the unification of**
  **biology**, *Nat. Genet.*, 25, 25-29.

- The Gene Ontology Consortium (2021) **The Gene Ontology resource: enriching**
  **a GOld mine**, *Nucleic Acids Res.*, 49, D325-D334.

- Gillespie, M. et al. (2022) **The reactome pathway knowledgebase 2022**,
  *Nucleic Acids Res.*, 50, D687-D692.

- Giurgiu, M et al. (2019) **CORUM: the comprehensive resource of mammalian**
  **protein complexes-2019**, *Nucleic Acids Res.*, 47, D559-D563

- Licata, L. et al. (2012) **MINT, the molecular interaction database: 2012**
  **update**, *Nucleic Acids Res.*, 40, D857-D861.

- Orchard, S. et al. (2014) **The MIntAct project-IntAct as a common curation**
  **platform for 11 molecular interaction databases**, *Nucleic Acids Res.*, 42,
  D358-D363.

- Oughtred, R. et al. (2018) **The BioGRID database: A comprehensive**
  **biomedical resource of curated protein, genetic, and chemical**
  **interactions**, *Protein Sci.*, 30, 187-200.

- Szklarczyk, D. et al. (2019) **STRING v11: protein-protein association**
  **networks with increased coverage, supporting functional discovery in**
  **genome-wide experimental datasets**, *Nucleic Acids Res.*, 47, D607-D613.

- The UniProt Consortium (2021) **UniProt: the universal protein knowledgebase**
  **in 2021**, *Nucleic Acids Res.*, 49, D480-D489.

---

The following applications are targeted.

- Shannon, P. et al. (2003) **Cytoscape: a software environment for integrated**
  **models of biomolecular interaction networks**, *Genome Res.*, 13, 2498-2504.

---

The following external libraries are utilized.

- Hagberg, A. A. et al. (2008) **Exploring network structure, dynamics, and**
  **function using NetworkX**, *Proceedings of the 7th Python in Science*
  *Conference*, 11-15.

- McKinney, W. (2010) **Data Structures for Statistical Computing in Python**,
  *Proceedings of the 9th Python in Science Conference*, 56-61.

- Virtanen, P. et al. (2020)  **SciPy 1.0: Fundamental Algorithms for**
  **Scientific Computing in Python**, *Nat. Methods*, 17, 261-272.

---

Development was inspired by previous work combining the following applications.

- Maere, S. et al. (2005) ***BiNGO*: a Cytoscape plugin to assess**
  **overrepresentation of Gene Ontology categories in Biological Networks**,
  *Bioinformatics*, 21, 3448-3449.

- Morris, J. H. et al. (2011) ***clusterMaker*: a multi-algorithm clustering**
  **plugin for Cytoscape**, *BMC Bioinform.*, 12.

- Su, G. et al. (2010) **GLay: community structure analysis of biological**
  **networks**, *Bioinformatics*, 26, 3135-3137.

---

References for implemented algorithms are supplied in the corresponding source
code.
