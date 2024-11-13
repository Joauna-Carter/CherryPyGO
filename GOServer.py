import cherrypy as cp
from Ontology import Ontology

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

if __name__ == '__main__':
    cp.quickstart(GOServer())
