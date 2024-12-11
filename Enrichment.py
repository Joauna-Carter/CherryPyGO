class Enrichment:
    def __init__(self, _term=None, _unadjPval=None, _pval=None, _genes=None):
        self.term = _term  # The OntoTerm
        self.unadjPval = _unadjPval  # Unadjusted p-value
        self.pval = _pval  # Adjusted p-value (FDR)
        self.genes = _genes if _genes is not None else set()  # Overlapping genes
        
    def toJSON(self):
        return {
            "term": self.term.toJSON(),
            "unadjusted_pvalue": self.unadjPval,
            "adjusted_pvalue": self.pval,
            "direct_annotations": len(self.term.directAnnos()),
            "total_annotations": len(self.term.allAnnos()),
            "overlapping_genes": [g.toJSON() for g in self.genes]
        }
