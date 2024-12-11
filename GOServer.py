import cherrypy as cp
from Ontology import Ontology
import cherrypy_cors
import json
from collections import defaultdict

cherrypy_cors.install()

#Update the CherryPy configuration
cp.config.update('server.conf')

class GOServer(object):

    def __init__(self):
        #Fields
        self.go = Ontology("go.obo","gene_association.sgd")

    @cp.expose
    def index(self, org=''):
        return "GO Server for organism " + org + " has " + str(len(self.go.roots)) + " roots with " + str(len(self.go.genes)) + " genes"

    @cp.expose
    @cp.tools.json_out()
    def gene(self, org='', uid=''):
        if(uid not in self.go.genes):
            return "Gene " + uid + " not found"
        else:
            g = self.go.genes[uid]
            msg =  {}
            msg["uid"] = g.uid
            msg["symbol"] = g.symbol
            msg["name"] = g.name
            msg["aliases"] = list(g.aliases)
            msg["terms"] = [x.toJSON() for x in g.terms()]
            return msg

    @cp.expose
    @cp.tools.json_out()
    def term(self, org='', uid=''):
        uid = uid.upper()
        if uid[0:3] != "GO:":
            uid = "GO:" + uid

        if(uid not in self.go.terms):
            return "Term " + uid + " not found"
        else:
            t = self.go.terms[uid]
            msg = {}
            msg["uid"] = t.uid
            msg["name"] = t.name
            msg["defn"] = t.defn
            msg["all"] = [x.toJSON() for x in list(t.allAnnos())]
            return msg

    @cp.expose
    @cp.tools.json_out()
    def allTerms(self, org=''):
        terms = []
        for termid in self.go.terms.keys():
            t = self.go.terms[termid]
            if len(t.allAnnos()) > 0:
                term = {}
                term["uid"] = t.uid
                term["name"] = t.name
                term["defn"] = t.defn
                term["all"] = [x.toJSON() for x in list(t.allAnnos())]
                terms.append(term)
        return terms

    @cp.expose
    @cp.tools.json_out()
    def tm(self, org='', uid=''):
        uid = uid.upper()
        if uid[0:3] != "GO:":
            uid = "GO:" + uid

        if (uid not in self.go.terms):
            return "Term " + uid + " not found"
        else:
            t = self.go.terms[uid]
            msg = t.toTreeMapJSON()
            return msg

    def _resolve_gene_names(self, query_list):
        """Resolve query strings to official gene symbols, handling aliases and duplicates.
        
        Returns:
            tuple containing:
            - dict mapping resolved queries to Gene objects
            - dict of name resolution details
        """
        resolution_results = {
            "direct_matches": [],    # Queries matching official symbols
            "alias_matches": [],     # Queries matching unique aliases
            "ambiguous_queries": [], # Queries matching multiple genes via aliases
            "unmatched_queries": []  # Queries with no matches
        }
        
        resolved_genes = {}  # Maps query string to Gene object
        
        # First pass: check for exact matches with official symbols
        for query in query_list:
            query = str(query).strip()
            if query in self.go.genes:
                resolved_genes[query] = self.go.genes[query]
                resolution_results["direct_matches"].append(query)
            else:
                # Build reverse alias lookup for unmatched queries
                matching_genes = []
                for gene in self.go.genes.values():
                    if query in gene.aliases:
                        matching_genes.append(gene)
                
                if len(matching_genes) == 1:
                    # Unique alias match
                    gene = matching_genes[0]
                    resolved_genes[query] = gene
                    resolution_results["alias_matches"].append({
                        "query": query,
                        "matched_to": gene.symbol
                    })
                elif len(matching_genes) > 1:
                    # Ambiguous alias
                    resolution_results["ambiguous_queries"].append({
                        "query": query,
                        "possible_matches": [g.symbol for g in matching_genes]
                    })
                else:
                    # No matches found
                    resolution_results["unmatched_queries"].append(query)
        
        return resolved_genes, resolution_results

    def _process_enrichment(self, org, genes, threshold):
        """Internal method to process enrichment analysis."""
        try:
            # Convert threshold to float
            threshold = float(threshold)
            
            # Ensure genes is a list
            if isinstance(genes, str):
                query_list = [g.strip() for g in genes.split(',') if g.strip()]
            else:
                query_list = [str(g).strip() for g in genes if str(g).strip()]
            
            if not query_list:
                return {"error": "No valid gene queries provided"}
            
            # Resolve gene names
            resolved_genes, resolution_results = self._resolve_gene_names(query_list)
            
            if not resolved_genes:
                return {
                    "error": "No valid genes found after name resolution",
                    "name_resolution": resolution_results
                }
            
            # Perform enrichment analysis with resolved genes
            results = self.go.enrichByGeneNames(list(resolved_genes.values()))
            
            # Filter by adjusted p-value threshold
            filtered_results = [r.toJSON() for r in results if r.pval <= threshold]
            
            return {
                "name_resolution": resolution_results,
                "query_summary": {
                    "input_count": len(query_list),
                    "resolved_count": len(resolved_genes),
                    "excluded_count": len(resolution_results["ambiguous_queries"]) + 
                                   len(resolution_results["unmatched_queries"])
                },
                "threshold": threshold,
                "total_results": len(results),
                "filtered_results": len(filtered_results),
                "enriched_terms": filtered_results
            }
            
        except ValueError:
            return {"error": "Invalid threshold value. Must be a number between 0 and 1"}
        except Exception as e:
            return {"error": str(e)}

    @cp.expose
    @cp.tools.json_out()
    def enrich(self, org='', genes='', threshold=0.05, **kwargs):
        """Perform enrichment analysis on a list of genes.
        Handles both GET and POST requests.
        
        For GET requests:
            genes: Comma-separated list of gene names/symbols/aliases
            threshold: P-value threshold for filtering results (default 0.05)
            
        For POST requests:
            Accepts JSON body with format:
            {
                "genes": ["gene1", "gene2", ...],  # List of gene names/symbols/aliases
                "threshold": 0.05
            }
        """
        # Check if this is a POST request with JSON data
        if cp.request.method == "POST":
            try:
                data = json.loads(cp.request.body.read().decode('utf-8'))
                genes = data.get('genes', [])
                threshold = data.get('threshold', 0.05)
            except json.JSONDecodeError:
                return {"error": "Invalid JSON in request body"}
            except Exception as e:
                return {"error": f"Error processing POST request: {str(e)}"}
        
        return self._process_enrichment(org, genes, threshold)

if __name__ == '__main__':
    print("Running main")
    cherrypy_cors.install()
    
    # Update config to handle POST requests with CORS
    config = {
        '/': {
            'cors.expose.on': True,
            'tools.response_headers.on': True,
            'tools.response_headers.headers': [
                ('Content-Type', 'application/json'),
                ('Access-Control-Allow-Origin', '*'),
                ('Access-Control-Allow-Methods', 'GET, POST, OPTIONS'),
                ('Access-Control-Allow-Headers', 'Content-Type'),
            ],
        }
    }
    cp.quickstart(GOServer(), config=config)
