# DIANA

DIANA is a command line application for **D**ata **I**ntegration **A**nd
**N**etwork-based **A**nalysis for post-translational modification mass
spectrometry data.

Incorporating protein-protein interactions from BioGRID, CORUM, IntAct, MINT,
Reactome, and STRING, DIANA assembles and analyses protein-protein interaction
networks from mass spectrometry data, providing integrative contextualization of
differential post-translational modification captured.

The enrichment of CORUM protein complexes, Gene Ontology terms, and Reactome
pathways by a protein-protein interaction network and its densely interacting
communities can be assessed as well as the distribution of mass spectrometric
measurements across communities.

Networks of Gene Ontology terms or Reactome pathways, reporting enrichment of
each term or pathway by a protein-protein interaction network, can be generated
along customized Cytoscape style specifications for each type of network.

## Setup
External dependencies can be installed using pip by running the following
command:

```
pip3 install -r diana/requirements.txt
```
External dependencies consist of NetworkX, pandas, and SciPy.

DIANA is currently developed using Python 3.10.4, Ubuntu 22.04 and Cytoscape
3.9.1.

## Command Line Interface

Instructions can be displayed by running the following command:

```
python3 diana/diana.py --help
```

Configuration is detailed below. Configuration files referring to published data
sets referenced below are included in this repository.

## Configuration

A configuration file specifies an array of workflows executed sequentially.
Multiple configuration files are processed concurrently.

---

Input genes or proteins are read from tabular input files.

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
from a cell entry in the table, possibly removing additional components of the
entry. The default setting is `"^(.+?)$"`, corresponding to the entire entry.

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
The spreadsheet from a file to extract input from. The default setting is `1`
corresponding to the first spreadsheet of the file.

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
The line number of the header, allowing to skip preceding lines. The default
setting is `1`, corresponding to the first line of the sheet.

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
A list of protein-protein interaction networks as exported from workflows. The
union of input protein-protein interaction networks with input proteins and
genes is used.

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
file for reference in the specification of further analyses. The default setting
is `0`.

```json
[
    {
      "proteins": [
        {
          "post-translational modification": "PTM"
        }
      ]
    }
]
```
The identifier for the type of post-translational modification associate with
measurements from the input file for reference in the specification of further
analyses. The default setting is `"PTM"`. Currently, up to two types of
post-translational modification per time of measurement are supported by
Cytoscape styles.

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
which they are ordered. If an entry contains fewer positions than measurements,
missing modification sites are substituted by 0. If an entry contains more
positions than measurements, only as many leading entries as there are
measurements are used.

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
A regular expression used to extract matching modification sites from an
entry in the table, allowing to remove additions to the site number. The default
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
A regular expression used to extract matching replicate measurements from an
entry in the table, allowing to remove additions to the measurement. The default
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
An inclusive threshold on the number of replicates required to consider a
measurement. The default setting is `1`.

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
largest absolute measurement. The default setting is `5`.

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
The function used to combine replicates into a single measurement. The function
is applied to ratios, not their binary logarithm. Cytoscape styles refer to this
combination. The default setting is `"mean"`. Available settings are `"mean"`,
`"median"`, `"max"`, `"maxabs"`, `"min"`, `"minabs"`, `"sum"`, and `"sumabs"`.

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
The base of the logarithm that measured measurements are reported as. By
default, ratios are assumed and no conversion is performed, corresponding to
`null`. Available settings are `null`, `2` and `10`.

---

The interface to sources of protein-protein interactions for the protein-protein
interaction network. The protein-protein interaction network is exported only if
any source is queried. Database-specific requirements can be defined, where
each must be satisfied for an interaction to be incorporated.

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
The function used to combine distinct modification sites into protein-specific
measurements. The function is applied to ratios, not their binary logarithm. The
default setting is `"maxabs"`, corresponding to the largest absolute
measurement. Available settings are `"mean"`, `"median"`, `"max"`, `"maxabs"`,
`"min"`, `"minabs"`, `"sum"` and `"sumabs"`.

```json
[
    {
      "Cytoscape": {
        "bar chart": {
          "replicate combination": "mean"
        },
        "node color": {
          "replicate combination": "mean"
        }
      }
    }
]
```
The function used to combine distinct replicates into modification site-specific
measurements. The function is applied to ratios, not their binary logarithm. The
default setting is `"mean"`, corresponding to the mean of replicates. Available
settings are `"mean"`, `"median"`, `"max"`, `"maxabs"`, `"min"`, `"minabs"`,
`"sum"` and `"sumabs"`.

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
The conversion that a measurement range refers to. It defaults to the binary
logarithm of a measurement for `null` but can be set to `"ratio"`,
`"percentile"`, `"quantile"` or `"standard score"`, computed with respect to the
distribution of a particular modification at a particular time of measurement
across the protein-protein interaction network.

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
The range of the bar charts reporting measurements. The adaptive default setting
is `[-1.0, 1.0]` if `"conversion"` is not set, `[0.5, 2.0]` if `"conversion"` is
set to `"ratio"`, `[25.0, 75.0]` if `"conversion"` is set to `"percentile"`,
`[0.25, 0.75]` if `"conversion"` is set to `"quantile"`, and `[-1.0, 1.0]` if
`"conversion"` is set to `"standard score"`.

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
exceeded or not. The adaptive default setting is `[-1.0, 1.0]` if `"conversion"`
is not set, `[0.5, 2.0]` if `"conversion"` is set to `"ratio"`, `[25.0, 75.0]`
if `"conversion"` is set to `"percentile"`, `[0.25, 0.75]` if `"conversion"` is
set to `"quantile"`, and `[-1.0, 1.0]` if `"conversion"` is set to
`"standard score"`.

```json
[
    {
      "Cytoscape": {
        "edge transparency": null,
      }
    }
]
```
The function used to combine edge confidence scores from in IntAct, MINT and
STRING, and, lacking a corresponding score, 1.0 for all interactions from
BioGRID, CORUM and Reactome. The combined score is reflected by edge
transparency in Cytoscape. By default, `null`, any edge receives a score of
1.0 and edges are not transparent. Available settings are `null`, `"mean"`,
`"median"`, `"max"`, `"min"`, `"sum"`,and `"number"`, the number of queried
databases supporting the protein-protein interaction.

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
The statistical test to assess enrichment of each complex, term or pathway by
the submitted proteins. The default setting is `"hypergeometric"`. Available
settings are `"binomial"` and `"hypergeometric"`.

```json
[
    {
      "CORUM enrichment": {
        "increase": true
      },
      "Gene Ontology enrichment": {
        "increase": true
      },
      "Reactome enrichment": {
        "increase": true
      }
    }
]
```
If `true`, assess overrepresentation, otherwise underrepresentation. The default
setting is `true`.

```json
[
    {
      "CORUM enrichment": {
        "correction": "Benjamini-Yekutieli"
      },
      "Gene Ontology enrichment": {
        "correction": "Benjamini-Yekutieli"
      },
      "Reactome enrichment": {
        "correction": "Benjamini-Yekutieli"
      }
    }
]
```
The procedure to correct p-values for multiple testing. The default setting is
`"Benjamini-Yekutieli"`. Available settings are `"Benjamini-Hochberg"`,
`"Benjamini-Yekutieli"`, `"Holm"` and `"Hommel"`.

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
The function used to combine distinct modification sites into protein-specific
measurements. The function is applied to ratios, not their binary logarithm.
The default setting is `"maxabs"`, corresponding to the largest absolute
measurement. Available settings are `"mean"`, `"median"`, `"max"`, `"maxabs"`,
`"min"`, `"minabs"`, `"sum"` and `"sumabs"`.

```json
[
    {
      "CORUM enrichment": {
        "subsets": [
          {
            "replicate combination": "mean"
          }
        ]
      },
      "Gene Ontology enrichment": {
        "subsets": [
          {
            "replicate combination": "mean"
          }
        ]
      },
      "Reactome enrichment": {
        "subsets": [
          {
            "replicate combination": "mean"
          }
        ]
      }
    }
]
```
The function used to combine distinct replicates into modification site-specific
measurements. The function is applied to ratios, not their binary logarithm.
The default setting is `"mean"`, corresponding to the mean of replicates.
Available settings are `"mean"`, `"median"`, `"max"`, `"maxabs"`, `"min"`,
`"minabs"`, `"sum"` and `"sumabs"`.

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
The conversion that a measurement range refers to. It defaults to the binary
logarithm of a measurement for `null` but can be set to `"ratio"`,
`"percentile"`, `"quantile"` or `"standard score"`, computed with respect to the
distribution of a particular modification at a particular time of measurement
across the protein-protein interaction network.

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
exceeded or not. The adaptive default setting is `[-1.0, 1.0]` if `"conversion"`
is not set, `[0.5, 2.0]` if `"conversion"` is set to `"ratio"`, `[25.0, 75.0]`
if `"conversion"` is set to `"percentile"`, `[0.25, 0.75]` if `"conversion"` is
set to `"quantile"`, and `[-1.0, 1.0]` if `"conversion"` is set to
`"standard score"`.

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
          "cellular component",
          "molecular function",
          "biological process"
        ]
      }
    }
]
```
The Gene Ontology namespaces to consider. The default setting is
`["cellular component", "molecular function" "biological process"]`.

---

Networks of Gene Ontology terms or Reactome pathways can be exported. Both
report the enrichment of each term or pathway by proteins represented in the
protein-protein interaction network with respect to the annotation specific to
an organism of interest. Edges of these networks are the respective hierarchical
relations of entities in these databases.

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
The function used to combine distinct modification sites into protein-specific
measurements. The function is applied to ratios, not their binary logarithm. The
default setting is `"maxabs"`, corresponding to the largest absolute
measurement. Available settings are `"mean"`, `"median"`, `"max"`, `"maxabs"`,
`"min"`, `"minabs"`, `"sum"` and `"sumabs"`.

```json
[
    {
      "Gene Ontology network": {
        "subsets": [
          {
            "replicate combination": "mean"
          }
        ]
      },
      "Reactome network": {
        "subsets": [
          {
            "replicate combination": "mean"
          }
        ]
      }
    }
]
```
The function used to combine distinct replicates into modification site-specific
measurements. The function is applied to ratios, not their binary logarithm. The
default setting is `"mean"`, corresponding to the mean of replicates. Available
settings are `"mean"`, `"median"`, `"max"`, `"maxabs"`, `"min"`, `"minabs"`,
`"sum"` and `"sumabs"`.

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
The conversion that a measurement range refers to. It defaults to the binary
logarithm of a measurement for `null` but can be set to `"ratio"`,
`"percentile"`, `"quantile"` or `"standard score"`, computed with respect to the
distribution of a particular modification at a particular time of measurement
across the protein-protein interaction network.

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
exceeded or not. The adaptive default setting is `[-1.0, 1.0]` if `"conversion"`
is not set, `[0.5, 2.0]` if `"conversion"` is set to `"ratio"`, `[25.0, 75.0]`
if `"conversion"` is set to `"percentile"`, `[0.25, 0.75]` if `"conversion"` is
set to `"quantile"`, and `[-1.0, 1.0]` if `"conversion"` is set to
`"standard score"`.

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
The statistical test to assess enrichment of each term or pathway by the
submitted proteins. The default setting is `"hypergeometric"`. Available
settings are `"binomial"` and `"hypergeometric"`.

```json
[
    {
      "Gene Ontology network": {
        "increase": true
      },
      "Reactome network": {
        "increase": true
      }
    }
]
```
If `true`, assess overrepresentation, otherwise underrepresentation. The default
setting is `true`.

```json
[
    {
      "Gene Ontology network": {
        "correction": "Benjamini-Yekutieli"
      },
      "Reactome network": {
        "correction": "Benjamini-Yekutieli"
      }
    }
]
```
The procedure to correct p-values for multiple testing. The default setting is
`"Benjamini-Yekutieli"`. Available settings are `"Benjamini-Hochberg"`,
`"Benjamini-Yekutieli"`, `"Holm"` and `"Hommel"`.

```json
[
    {
      "Gene Ontology network": {
        "annotation": false
      },
      "Reactome network": {
        "annotation": false
      }
    }
]
```
If `true`, compute enrichment with respect to the entire annotation, specific to
the organism of interest, otherwise with respect to proteins represented in the
protein-protein interaction network to assess enrichment within specified
subsets of them. The default setting is `false`.

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
          "cellular component",
          "molecular function",
          "biological process"
        ]
      }
    }
]
```
The Gene Ontology namespaces to consider. The default setting is
`["cellular component", "molecular function" "biological process"]`.

---

Communities of the protein-protein interaction network can be extracted using
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
The function used to combine edge confidence scores from in IntAct, MINT and
STRING, and, lacking of corresponding score, 1.0 for all interactions from
BioGRID, CORUM and Reactome. The combined score is used as edge weight in
community detection. By default, `null`, any edge receives a score of 1.0,
corresponding to an unweighted network. Available settings are `null`, `"mean"`,
`"median"`, `"max"`, `"min"`, `"sum"`, and `"number"`, the number of queried
databases supporting the protein-protein interaction.

```json
[
    {
      "community detection": {
        "community size": null
      }
    }
]
```
An upper bound on the number of proteins per community. Modules are iteratively
subdivided until this threshold is met. The adaptive default setting is the
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
The function used to combine community sizes in terms of nodes into the value
decisive to meeting the community size threshold. The default setting is
`"mean"`. Available settings are `"mean"`, `"median"`, `"max"`, and `"min"`.

---

CORUM protein complex, Gene Ontology term and Reactome pathway enrichment by
separate communities can be assessed. The proteins considered can be
restricted, based on ranges of associated mass spectrometric measurements,
either by a union or intersection of specified subsets of proteins from the
protein-protein interaction network.

To assess the distribution of mass spectrometry measurements across communities,
these can be categorized to measure the enrichment of proteins exceeding a
specified absolute or relative threshold.

Alternatively, the distribution of measurements within separate communities can
be compared with the remaining network.

A community is exported if is significant according to any of the specified
tests.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "increase": true
        },
        "Gene Ontology enrichment": {
          "increase": true
        },
        "Reactome enrichment": {
          "increase": true
        },
        "measurement enrichment": {
          "increase": true
        },
        "measurement location": {
          "increase": true
        }
      }
    }
]
```
If `true`, assess overrepresentation, otherwise underrepresentation concerning
enrichment and relative increase, otherwise relative decrease of measurements in
a community. The default setting is `true`.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "correction": "Benjamini-Yekutieli"
        },
        "Gene Ontology enrichment": {
          "correction": "Benjamini-Yekutieli"
        },
        "Reactome enrichment": {
          "correction": "Benjamini-Yekutieli"
        },
        "measurement enrichment": {
          "correction": "Benjamini-Yekutieli"
        },
        "measurement location": {
          "correction": "Benjamini-Yekutieli"
        }
      }
    }
]
```
The procedure to correct p-values for multiple testing. The default setting is
`"Benjamini-Yekutieli"`. Available settings are `"Benjamini-Hochberg"`,
`"Benjamini-Yekutieli"`, `"Holm"` and `"Hommel"`.

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
The statistical test to assess enrichment of each complex, term or pathway by
the each of the communities. The default setting is `"hypergeometric"`.
Available settings are `"binomial"` and `"hypergeometric"`.

```json
[
    {
      "community detection": {
        "measurement location": {
          "test": "Mann-Whitney-Wilcoxon"
        }
      }
    }
]
```
The statistical test to compare locations of modification- and time-specific
distributions of measurements across each community with the distribution of the
remaining protein-protein interaction network. The default setting is
`"Mann-Whitney-Wilcoxon"`. Available settings are `"Mann-Whitney-Wilcoxon"` and
`"Welch"`.

```json
[
    {
      "community detection": {
        "measurement location": {
          "absolute": true
        }
      }
    }
]
```
If `true` test the absolute values of the measurements, otherwise the
measurements. The default value is `true`.

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
the organism of interest, otherwise with respect to proteins represented in the
protein-protein interaction network. The default setting is `false`.

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
The function used to combine distinct modification sites into protein-specific
measurements. The function is applied to ratios, not their binary logarithm. The
default setting is `"maxabs"`, corresponding to the largest absolute
measurement. Available settings are `"mean"`, `"median"`, `"max"`, `"maxabs"`,
`"min"`, `"minabs"`, `"sum"` and `"sumabs"`.

```json
[
    {
      "community detection": {
        "CORUM enrichment": {
          "subsets": [
            {
              "replicate combination": "mean"
            }
          ]
        },
        "Gene Ontology enrichment": {
          "subsets": [
            {
              "replicate combination": "mean"
            }
          ]
        },
        "Reactome enrichment": {
          "subsets": [
            {
              "replicate combination": "mean"
            }
          ]
        }
      }
    }
]
```
The function used to combine distinct replicates into modification site-specific
measurements. The function is applied to ratios, not their binary logarithm. The
default setting is `"mean"`, corresponding to the mean of replicates. Available
settings are `"mean"`, `"median"`, `"max"`, `"maxabs"`, `"min"`, `"minabs"`,
`"sum"` and `"sumabs"`.

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
The conversion that a measurement range refers to. It defaults to the binary
logarithm of a measurement for `null` but can be set to `"ratio"`,
`"percentile"`, `"quantile"` or `"standard score"`, computed with respect to the
distribution of a particular modification at a particular time of measurement
across the protein-protein interaction network. Conversions refer to separate
communities.

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
exceeded or not. The adaptive default setting is `[-1.0, 1.0]` if `"conversion"`
is not set, `[0.5, 2.0]` if `"conversion"` is set to `"ratio"`, `[25.0, 75.0]`
if `"conversion"` is set to `"percentile"`, `[0.25, 0.75]` if `"conversion"` is
set to `"quantile"`, and `[-1.0, 1.0]` if `"conversion"` is set to
`"standard score"`.

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
              "cellular component",
              "molecular function",
              "biological process"
            ]
        }
      }
    }
]
```
The Gene Ontology namespaces to consider. The default setting is
`["cellular component", "molecular function" "biological process"]`.

```json
[
    {
      "community detection": {
        "measurement enrichment": {
            "site combination": "maxabs"
        },
        "measurement location": {
            "site combination": "maxabs"
        }
      }
    }
]
```
The function used to combine distinct modification sites into protein-specific
measurements. The default setting is `"maxabs"`, corresponding to the largest
absolute value. Available settings are `"mean"`, `"median"`, `"max"`,
`"maxabs"`, `"min"`, `"minabs"`, `"sum"`, `"sumabs"` and `null` to consider
modification sites separately.

```json
[
    {
      "community detection": {
        "measurement enrichment": {
            "replicate combination": "mean"
        },
        "measurement location": {
            "replicate combination": "mean"
        }
      }
    }
]
```
The function used to combine distinct replicates into modification site-specific
measurements. The default setting is `"mean"`, corresponding to the mean of
replicates. Available settings are `"mean"`, `"median"`, `"max"`, `"maxabs"`,
`"min"`, `"minabs"`, `"sum"`, `"sumabs"` and `null` to consider replicates
separately.

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
The conversion that a measurement range refers to. It defaults to the binary
logarithm of a measurement for `null` but can be set to `"ratio"`,
`"percentile"`, `"quantile"` or `"standard score"`, computed with respect to the
distribution of a particular modification at a particular time of measurement
across each community of the protein-protein interaction network.

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
The range of measurements to group proteins by whether the range is exceeded or
not. The adaptive default setting is `[-1.0, 1.0]` if `"conversion"` is not set,
`[0.5, 2.0]` if `"conversion"` is set to `"ratio"`, `[25.0, 75.0]` if
`"conversion"` is set to `"percentile"`, `[0.25, 0.75]` if `"conversion"` is set
to `"quantile"`, and `[-1.0, 1.0]` if `"conversion"` is set to
`"standard score"`.

## Protein-protein interaction network

Annotations of proteins contain the following information:

```xml
<node id="Q12933">
  <data key="30 U 1 1">3.765004</data>
  <data key="30 U 1 2">4.2322</data>
  <data key="30 U 2 1">3.926948</data>
  <data key="30 U 2 2">4.299611</data>
  <data key="30 U 2 3">3.715125</data>
  <data key="30 U 3 1">2.239367</data>
  <data key="30 U 3 2">3.353323</data>
  <data key="30 U 4 1">4.508999</data>
  <data key="30 U 4 2">4.576885</data>
  <data key="120 U 1 1">2.377096</data>
  <data key="120 U 1 2">2.890972</data>
  <data key="120 U 2 1">2.962512</data>
  <data key="120 U 2 2">3.0054</data>
  <data key="gene">TRAF2</data>
  <data key="protein">TNF receptor-associated factor 2</data>
  <data key="post-translational modification 30">U</data>
  <data key="post-translational modification 120">U</data>
  <data key="30 U 1">4.017431773281321</data>
  <data key="30 U 2">4.001081644434266</data>
  <data key="30 U 3">2.9012913530601145</data>
  <data key="30 U 4">4.543341260044583</data>
  <data key="measurement 30">UP</data>
  <data key="120 U 1">2.6567938583308726</data>
  <data key="120 U 2">2.9841153643117204</data>
  <data key="measurement 120">UP</data>
</node>
```
Proteins are represented by their primary UniProt accession. `"gene"` and
`"protein"` refer to the gene and protein names listed in UniProt, respectively.

Extracted mass spectrometric measurements, expressed as binary logarithms of the
corresponding ratios, are represented by keys consisting of four components.
The first component refers to the specified time of measurement, the second to
the specified type of post-translational modification, the third to the relative
position of the corresponding modification site and the fourth to the replicate.
Entries consisting of three components refer to modification site-specific
combinations of replicates.

`"post-translational modification 30"` and
`"post-translational modification 120"` each refer to the specified types of
post-translational modification the protein is associated with while
`"measurement 30"` and `"measurement 120"` each refer to a categorization of the
protein-specific combination of modification sites for particular specified
times of measurement.

Annotations of protein-protein interactions contain the following information:

```xml
<edge source="Q92844" target="Q12933">
  <data key="BioGRID">1.0</data>
  <data key="CORUM">1.0</data>
  <data key="IntAct">0.81</data>
  <data key="Reactome">1.0</data>
  <data key="STRING">0.991</data>
  <data key="score">1.0</data>
</edge>
```
A protein-protein interaction between is associated with database-specific
confidence scores. `"score"` refers to the combination of confidence scores.

## Gene Ontology network

Annotations of Gene Ontology terms contain the following information:

```xml
<node id="GO:0043122">
  <data key="term">regulation of I-kappaB kinase/NF-kappaB signaling</data>
  <data key="namespace">biological process</data>
  <data key="p-value">0.07496498697418053</data>
  <data key="number of proteins">5</data>
  <data key="proteins">P04792 P55085 Q12933 Q13501 Q9Y4K3</data>
</node>
```
Terms are represented by their Gene Ontology ID. `"term"` refers to the term and
`"namespace"` to its namespace. `"p-value"` refers to the corrected p-value of
the test for enrichment of the term by the submitted proteins.
`"number of proteins"` refers to the number of submitted proteins annotated with
the term and `"proteins"` refers to the space separated accessions of these
proteins.

Relationships of Gene Ontology terms contain no additional information besides
the related terms.

## Reactome network

Annotations of Reactome pathways contain the following information:

```xml
<node id="R-HSA-9758274">
  <data key="pathway">Regulation of NF-kappa B signaling</data>
  <data key="p-value">1.0</data>
  <data key="number of proteins">3</data>
  <data key="proteins">O75113 Q12933 Q9Y4K3</data>
</node>
```
Pathways are represented by their Reactome pathway stable identifier.
`"pathway"` refers to the pathway name. `"p-value"` refers to the corrected
p-value of the test for enrichment of the pathway by the submitted proteins.
`"number of proteins"` refers to the number of submitted proteins annotated with
the pathway and `"proteins"` refers to the space separated accessions of these
proteins.

Relationships of Reactome pathways contain no additional information besides the
related pathways.

## References

The configuration files in this repository refer to data sets supplemented with
the following publications:

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

DIANA accesses the following resources:

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

DIANA utilizes the following external libraries:

- Hagberg, A. A. et al. (2008) **Exploring network structure, dynamics, and**
  **function using NetworkX**, *Proceedings of the 7th Python in Science*
  *Conference*, 11-15.

- McKinney, W. (2010) **Data Structures for Statistical Computing in Python**,
  *Proceedings of the 9th Python in Science Conference*, 56-61.

- Virtanen, P. et al. (2020)  **SciPy 1.0: Fundamental Algorithms for**
  **Scientific Computing in Python**, *Nat. Methods*, 17, 261-272.

---

Development of DIANA was inspired by previous work combining the following
applications:

- Maere, S. et al. (2005) ***BiNGO*: a Cytoscape plugin to assess**
  **overrepresentation of Gene Ontology categories in Biological Networks**,
  *Bioinformatics*, 21, 3448-3449.

- Morris, J. H. et al. (2011) ***clusterMaker*: a multi-algorithm clustering**
  **plugin for Cytoscape**, *BMC Bioinform.*, 12.

- Su, G. et al. (2010) **GLay: community structure analysis of biological**
  **networks**, *Bioinformatics*, 26, 3135-3137.

---

DIANA produces Cytoscape style specifications.

- Shannon, P. et al. (2003) **Cytoscape: a software environment for integrated**
  **models of biomolecular interaction networks**, *Genome Res.*, 13, 2498-2504.

---

References for implemented algorithms are listed in the corresponding source
code.

---

DIANA is developed by Jens Rieser and Lucas Fein in the Molecular Bioinformatics
group of Ina Koch at Goethe-University Frankfurt.