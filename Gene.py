from OntoTerm import OntoTerm
import json

class Gene:
    def __init__(self, _uid=None, _symb=None, _name=None, _alias=set()):
        self.uid = _uid
        self.symbol = _symb
        self.name = _name
        self.aliases = _alias
        self.direct = {}  # Code --> Set(OntoTerm)
        self.annos = {}  # Code --> Set(OntoTerm)
    def toJSON(self):
        return self.uid + "(" + self.symbol + ")"

    def terms(self, exclude=None):
        if exclude == None:
            exclude = []
        result = set()
        for code in self.annos:
            # If this is not an excluded evidence code
            if code not in exclude:
                result = result.union(self.annos[code])
        return result

    def directAnnotate(self, term, code):
        if code not in self.direct:
            self.direct[code] = set()
        self.direct[code].add(term)

    def addAnnos(self, term, code):
        if code not in self.annos:
            self.annos[code] = set()
        self.annos[code].add(term)
