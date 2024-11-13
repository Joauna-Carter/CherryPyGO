from OntoTerm import OntoTerm
from Gene import Gene

class Ontology:
    def __init__(self, oboFileName=None, annoFileName=None):
        #Fields
        self.terms = {} #Term.uid --> Term
        self.roots = [] #Root terms
        self.genes = {} #Gene.uid --> Gene

        #Load files if able
        if oboFileName is not None:
            self.loadOboFile(oboFileName)
            if annoFileName is not None:
                self.loadAnnoFile(annoFileName)


    def loadOboFile(self, oboFileName):
        fin = open(oboFileName,"r")
        #Parsing helper variables
        interm = False
        uid = ""
        name = ""
        defn = ""
        is_a = set()
        part_of = set()
        regulates = set()
        #Parse the OBO file format
        for line in fin:
            line = line.strip()
            if line == "":
                if interm: #If blank line and was reading in a term, add it
                    if uid not in self.terms:
                        self.terms[uid] = OntoTerm(uid,name,defn,is_a)
                    else:
                        term = self.terms[uid]
                        term.name = name
                        term.defn = defn
                        term.is_a = is_a
                        term.part_of = part_of
                        term.regulates = regulates
                    if len(self.terms[uid].parents()) == 0:
                        self.roots.append(self.terms[uid])
                interm = False
            elif line == "[Term]":
                #Reset everything for this term
                interm = True
                uid = ""
                name = ""
                defn = ""
                is_a = set()
                part_of = set()
                regulates = set()
            elif interm:
                if line[0:3] == "id:":
                    uid = line[4:]
                elif line[0:5] == "name:":
                    name = line[6:]
                elif line[0:4] == "def:":
                    defn = line[5:]
                elif line[0:12] == "is_obsolete:":
                    interm = False
                elif line[0:5] == "is_a:":
                    termid = line[6:16]
                    if termid not in self.terms:
                        self.terms[termid] = OntoTerm(termid)
                    is_a.add(self.terms[termid])
                elif line[0:21] == "relationship: part_of":
                    termid = line[22:32]
                    if termid not in self.terms:
                        self.terms[termid] = OntoTerm(termid)
                    part_of.add(self.terms[termid])
                elif line[0:23] == "relationship: regulates":
                    termid = line[24:34]
                    if termid not in self.terms:
                        self.terms[termid] = OntoTerm(termid)
                    regulates.add(self.terms[termid])
                elif line[0:34] == "relationship: positively_regulates" or line[0:34] == "relationship: negatively_regulates":
                    termid = line[35:45]
                    if termid not in self.terms:
                        self.terms[termid] = OntoTerm(termid)
                    regulates.add(self.terms[termid])
        fin.close()
        for term in self.terms.values():
            for parent in term.parents():
                parent.children.add(term)
                print(str(len(term.children)) + " ~~ " + str(len(parent.children)))
        print(len(self.terms["GO:0042867"].children))

    def loadAnnoFile(self, annoFileName):
        fin = open(annoFileName, "r")
        for line in fin:
            if line[0] != "!":
                line = line.strip()
                parts = line.split("\t")
                #Get gene instance (Create it if needed)
                if parts[2] not in self.genes:
                    self.genes[parts[2]] = Gene(parts[1], parts[2], parts[9], set(parts[10].split("|")))
                g = self.genes[parts[2]]
                #Annotate the gene to Terms
                if "NOT" not in parts[3]: #Eliminate negative annotations
                    if parts[4] not in self.terms:
                        print("ERROR: " + parts[4] + "not found in Ontology")
                    else:
                        term = self.terms[parts[4]]
                        term.annotate(g, parts[6])

    #Return map of term->pvalue for enrichment results
    def pvalsByGene(self, qGenes):
        termsToCheck = set()
        pvalResults = {}
        for gene in qGenes:
            termsToCheck = termsToCheck.union(gene.terms())
        for term in termsToCheck:
            pvalResults[term] = term.hyperGeomEnrich(qGenes, len(self.roots[0].allAnnos()))
        return pvalResults

    #Return map of term->pvalue for enrichment results
    def pvalsByGeneNames(self, qGeneNames):
        qGenes = set()
        for name in qGeneNames:
            if name in self.genes:
                qGenes.add(self.genes[name])
        if len(qGenes) > 0:
            return self.pvalsByGene(qGenes)
        else:
            return {}

    #Return list of enrichment records from enrichment search
    def enrichByGeneNames(self, qGeneNames):
        qGenes = set()
        for name in qGeneNames:
            if name in self.genes:
                qGenes.add(self.genes[name])
        if len(qGenes) == 0:
            return []
        else:
            termsToCheck = set()
            enrichResults = []
            for gene in qGenes:
                termsToCheck = termsToCheck.union(gene.terms())
            #for term in termsToCheck:
            #    enrichResults.append(term.enrichRecord(qGenes, /**TODO**/))


if __name__ == '__main__':
    go = Ontology("go.obo","gene_association.mgi")
    for root in go.roots:
        print("*" + root.name)
        print("**" + str(len(root.directAnnos())))
        print("**" + str(len(root.allAnnos())))
    print(len(go.genes))

    pvals = go.pvalsByGeneNames(["Notch1","Notch2","Notch3"])
    for term in pvals:
        print(term.uid + "\t" + str(pvals[term]) + "\t" + term.name)
