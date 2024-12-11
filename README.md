# CherryPyGO

A Python-based REST API for Gene Ontology (GO) enrichment analysis, built with CherryPy.

## Overview

This API provides access to Gene Ontology data and enrichment analysis capabilities. It handles gene name resolution (including aliases), performs hypergeometric enrichment tests, and applies FDR correction for multiple testing.

## Installation

1. Clone the repository
2. Install dependencies:
```bash
pip install cherrypy cherrypy_cors scipy
```
3. Place your GO data files:
   - `go.obo`: The Gene Ontology structure file
   - `gene_association.sgd`: Gene annotations file (SGD format)

## API Endpoints

### GET /
Returns basic statistics about the loaded GO data.

**Parameters:**
- `org`: Organism identifier (optional)

**Response:**
```json
"GO Server for organism [org] has [X] roots with [Y] genes"
```

### GET /gene
Retrieves information about a specific gene.

**Parameters:**
- `org`: Organism identifier (optional)
- `uid`: Gene identifier

**Response:**
```json
{
    "uid": "YBR020W",
    "symbol": "GAL1",
    "name": "Galactokinase",
    "aliases": ["GAL1", "GAL3"],
    "terms": [
        {
            "uid": "GO:0006012",
            "name": "galactose metabolic process",
            "defn": "..."
        }
    ]
}
```

### GET /term
Retrieves information about a specific GO term.

**Parameters:**
- `org`: Organism identifier (optional)
- `uid`: GO term ID (with or without 'GO:' prefix)

**Response:**
```json
{
    "uid": "GO:0006012",
    "name": "galactose metabolic process",
    "defn": "...",
    "all": ["YBR020W(GAL1)", "YBR019C(GAL10)"]
}
```

### GET /enrich
Performs GO enrichment analysis on a set of genes.

**Parameters:**
- `org`: Organism identifier (optional)
- `genes`: Comma-separated list of gene names/symbols/aliases
- `threshold`: P-value threshold for significance (default: 0.05)

**Response:**
```json
{
    "name_resolution": {
        "direct_matches": ["YBR020W", "GAL1"],
        "alias_matches": [
            {
                "query": "GAL10",
                "matched_to": "YBR019C"
            }
        ],
        "ambiguous_queries": [
            {
                "query": "ADH",
                "possible_matches": ["ADH1", "ADH2", "ADH3", "ADH4"]
            }
        ],
        "unmatched_queries": ["INVALID1"]
    },
    "query_summary": {
        "input_count": 5,
        "resolved_count": 3,
        "excluded_count": 2
    },
    "threshold": 0.05,
    "total_results": 150,
    "filtered_results": 12,
    "enriched_terms": [
        {
            "term": {
                "uid": "GO:0006012",
                "name": "galactose metabolic process",
                "defn": "..."
            },
            "unadjusted_pvalue": 0.0001,
            "adjusted_pvalue": 0.005,
            "direct_annotations": 5,
            "total_annotations": 15,
            "overlapping_genes": ["YBR020W(GAL1)", "YBR019C(GAL10)"]
        }
    ]
}
```

### POST /enrich
Same as GET /enrich but accepts larger gene lists via POST request.

**Request Body:**
```json
{
    "genes": ["GAL1", "GAL10", "ADH", "INVALID1"],
    "threshold": 0.05
}
```

**Response:** Same format as GET /enrich

## Gene Name Resolution

The API handles various ways that genes might be specified:
1. Official symbols (e.g., "GAL1")
2. Systematic names (e.g., "YBR020W")
3. Aliases/common names

The name resolution process:
1. First attempts exact matches to official symbols
2. For unmatched queries, searches through gene aliases
3. Handles ambiguous cases (e.g., one alias matching multiple genes)
4. Reports unmatched queries

## Statistical Methods

### Enrichment Analysis
- Uses hypergeometric test to calculate enrichment
- Applies Benjamini-Hochberg FDR correction for multiple testing
- Returns both raw and adjusted p-values
- Filters results based on user-specified significance threshold

### Annotation Propagation
GO term annotations are propagated up the ontology hierarchy following the true path rule:
- Direct annotations are those explicitly stated in the annotation file
- Total annotations include both direct and inherited annotations

## Error Handling

The API returns clear error messages for common issues:
- Invalid gene names
- Malformed requests
- Invalid threshold values
- Empty query lists

Example error response:
```json
{
    "error": "No valid genes found after name resolution",
    "name_resolution": {
        // Resolution details included for debugging
    }
}
```

## Running the Server

Start the server:
```bash
python GOServer.py
```

The server will start on port 8001 by default (configurable in server.conf).

## CORS Support

The API supports Cross-Origin Resource Sharing (CORS) and can be accessed from web applications on different domains.