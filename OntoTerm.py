import scipy.stats as stats

class OntoTerm:
    # uid = ""
    # name = ""
    # defn = ""
    # is_a = set() #Set of OntoTerms that are direct parents
    # part_of = set()
    # regulates = set()
    # children = set()

    def __init__(self, _uid="", _name="", _defn="", _is_a=set(), _part_of=set(), _regulates=set()):
        self.uid  = _uid
        self.name = _name
        self.defn = _defn
        self.is_a = _is_a
        self.part_of = _part_of
        self.regulates = _regulates
        self.direct = {} #Code --> Gene
        self.annos = {} #Code --> Gene
        self.children = set()

    def toJSON(self):
        return {"uid": self.uid,
                "name": self.name,
                "defn": self.defn }

    def toTreeMapPlotlyJSON(self, depth=2):
        labels = [self.uid]
        parents = [""]
        if len(self.parents()) > 0:
            parents = [list(self.parents())[0].uid]
        if depth > 0:
            for child in self.children:
                kidStuff = child.toTreeMapPlotlyJSON(depth - 1)
                for kidLabel in kidStuff["labels"]:
                    labels.append(kidLabel)
                for kidParent in kidStuff["parents"]:
                    parents.append(kidParent)
        #print(len(labels))
        #print(len(parents))
        #print(labels)
        return {
            'labels': labels,
            'parents': parents,
        }

    def toTreeMapJSON(self, depth=6, withAnno=True):
        data = { 'name': self.name, 'uid': self.uid }
        data['size'] = len(self.allAnnos())
        data['kids'] = []
        if depth > 0:
            for kid in self.children:
                if not withAnno or len(kid.allAnnos()) > 0:
                    data['kids'].append(kid.toTreeMapJSON(depth-1))
        return data

    def parents(self):
        return(self.is_a.union(self.part_of).union(self.regulates))

    def directAnnos(self):
        result = set()
        for code in self.direct:
            result = result.union(self.direct[code])
        return result

    def allAnnos(self):
        result = set()
        for code in self.annos:
            result = result.union(self.annos[code])
        return result

    def annotate(self, gene, code):
        #Add gene to term
        if code not in self.direct:
            self.direct[code] = set()
        self.direct[code].add(gene)
        #Add term to gene
        gene.directAnnotate(self,code)
        #Recursively annotate up
        self.addAnnos(gene, code)

    def addAnnos(self, gene, code):
        #Add gene to the annotations of this term
        if code not in self.annos:
            self.annos[code] = set()
        if gene not in self.annos[code]:
            self.annos[code].add(gene)
            #Add term to gene
            gene.addAnnos(self,code)
            #Recursively annotate the gene upward
            for parent in self.parents():
                parent.addAnnos(gene, code)

    def hyperGeomEnrich(self, genes, totalGenes):
        termGenes = self.allAnnos()
        hgd = stats.hypergeom(totalGenes, len(termGenes), len(genes))
        overlap = termGenes.intersection(genes)
        p = 0
        for x in range(len(overlap),len(genes)+1):
            p += hgd.pmf(x)
        return p


testOntolgoy = OntoTerm("GO:0000015",
                        "phosphopyruvate hydratase complex",
                        "A multimeric enzyme complex, usually a dimer or an octamer, that catalyzes the conversion of 2-phospho-D-glycerate to phosphoenolpyruvate and water.",
                        set(["GO:0044445","GO:1902494"]))
print(testOntolgoy)
