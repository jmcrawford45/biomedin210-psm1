import xml.etree.ElementTree as ET
import requests
import time
import nltk
import nltk.tag, nltk.data
from collections import Counter
import string
import os.path

BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/'
TOOL_IDENTIFIER = '&tool=hp-ontology&email=jared13@stanford.edu'
SEARCH_URL = BASE_URL + 'entrez/eutils/esearch.fcgi?db=pubmed&retmode=xml&retmax=500&term=hypersensitivity%+pneumonitis'
FETCH_URL = BASE_URL + 'entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id='
ANNOTATION_URL = 'https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/RESTful/tmTool.cgi/Chemical,Species/{}/BioC'
OCCUPATIONS_FILE = 'occupations.txt'
STOP_TERMS = ['DT', 'CC', ',', '.', ':', '!', 'CD', 'IN', 'PRP', 'PRP$', 'EX' , 'POS']
DEFAULT_TAGGER = nltk.data.load('taggers/maxent_treebank_pos_tagger/english.pickle')
PRINTABLE = frozenset(string.printable)
# Ignore humans as causative factors
SPECIES_STOP_WORDS = frozenset([
'human',
'patient',
'patients',
'woman',
'man',
'people',
'men',
'women'
])

#simple rate limiter algorithm found at
#https://stackoverflow.com/questions/667508/whats-a-good-rate-limiting-algorithm
def RateLimited(maxPerSecond):
    minInterval = 1.0 / float(maxPerSecond)
    def decorate(func):
        lastTimeCalled = [0.0]
        def rateLimitedFunction(*args,**kargs):
            elapsed = time.clock() - lastTimeCalled[0]
            leftToWait = minInterval - elapsed
            if leftToWait>0:
                time.sleep(leftToWait)
            ret = func(*args,**kargs)
            lastTimeCalled[0] = time.clock()
            return ret
        return rateLimitedFunction
    return decorate

class Article(object):

	def __init__(self, title, abstract, pmid):
		self.title = title
		self.abstract = abstract
		self.pmid = pmid
		self.text = title + '\n\n' + abstract

class OccupationResult(object):

	def __init__(self, name, pmid):
		self.names = set()
		self.name = name
		self.names.add(name)
		self.pmid = pmid
		self.causative_factors = set()

	def add_factor(self, factor):
		if factor.name.lower() not in SPECIES_STOP_WORDS:
			self.causative_factors.add(factor)

	def add_name(self, name):
		self.names.add(name)

	def __str__(self):
		return ('----------\nPMID: {}\n Occupations: {}\n Supsected Causative Factors\n [\n{}\n]\n---------'.format(
			self.pmid, self.names, '\n'.join([str(c) for c in self.causative_factors])))

class CausativeFactor(object):

	def __init__(self, cause_type, name, mesh):
		self.type = cause_type
		self.name = name
		self.mesh = mesh

	def __eq__(self, other):
		return self.type == other.type and self.name == other.name and self.mesh == other.mesh

	def __str__(self):
		return 'CausativeFactor type: {}, name: {}, mesh: {}'.format(
			self.type, self.name, self.mesh)

class HPExtractor(object):

	def __init__(self):
		occupations = self.get_raw_occupations()
		model = self.get_tagged_occupations(occupations)
		self.tagger = nltk.tag.UnigramTagger(model=model, backoff=DEFAULT_TAGGER)
		self.articles = list()
		self.occupations = list()

	def get_tagged_occupations(self, occupations):
		tokens = nltk.word_tokenize(occupations)
		tagged = nltk.pos_tag(tokens)
		tagged = [(word, 'OCC') for (word, tag) in tagged if tag not in STOP_TERMS and len(word) > 3]
		return dict(tagged)


	def get_raw_occupations(self):
		return open(OCCUPATIONS_FILE).read().lower()

	#Returns the root of an XML tree corresponding to the url parameter
	#1 RPS so we don't make anyone upset
	@RateLimited(1)
	def get_http(self, url, footer=True):
		if footer:
			url += TOOL_IDENTIFIER
		text = requests.get(url).text
		text = ''.join([c for c in text if c in PRINTABLE])
		return text

	def get_xml(self, url, footer=True):
		return ET.fromstring(self.get_http(url, footer))

	# Queries E-util search API for `retmax` most relevant articles
	def get_pmid_list(self):
		ids_list = list()
		root = self.get_xml(SEARCH_URL)
		for pmid in root.iter('Id'):
			ids_list.append(pmid.text)
		self.ids_list = ids_list

	def get_cached_article(self, article):
		fname = 'articles/{}.txt'.format(article)
		if os.path.isfile(fname):
			return ET.fromstring(open(fname).read())
		text = self.get_http(FETCH_URL+article)
		with open(fname, "w") as f:
			f.write(text)
		return ET.fromstring(text)	

	def fetch_articles(self):
		self.get_pmid_list()
		for article in self.ids_list:
			print 'Get article {}'.format(article)
			try:
				abstract = ''
				title = ''
				article_root = self.get_cached_article(article)
				for a in article_root.iter('AbstractText'):
					abstract = a.text
				for t in article_root.iter('ArticleTitle'):
					title = t.text
				self.articles.append(Article(abstract, title, article))
			except ET.ParseError as e:
				print e
				print FETCH_URL + article
				print 'Error parsing article with PMID={}'.format(article)
				continue
		print '{} articles loaded into extractor'.format(len(self.articles))

	def indentify_causes(self):
		articles = ','.join(self.occupations.keys())
		collection = self.get_xml(ANNOTATION_URL.format(articles), False)
		for article in collection.iter('document'):
			pmid = article.find('id').text
			for annotation in article.iter('annotation'):
				cause_type = ''
				mesh = ''
				name = annotation.find('text').text
				for infon in annotation.iter('infon'):
					if infon.attrib['key'] == 'type':
						cause_type = infon.text
					elif infon.attrib['key'] == 'MESH':
						mesh = infon.text
				factor = CausativeFactor(cause_type, name, mesh)
				self.occupations[pmid].add_factor(factor)

	def run(self):
		self.fetch_articles()
		for article in self.articles:
			tagged = self.tagger.tag(nltk.word_tokenize(article.text))
			#remove duplicates
			tagged = list(set(tagged))
			filtered = [OccupationResult(w, article.pmid) for (w,t) in tagged if t=='OCC']
			self.occupations += filtered
		occupation_dict = dict()
		for o in self.occupations:
			if o.pmid in occupation_dict:
				occupation_dict[o.pmid].add_name(o.name)
			else:
				occupation_dict[o.pmid] = o
		self.occupations = occupation_dict
		self.indentify_causes()
		for v in self.occupations.values():
			print str(v)


hp = HPExtractor()
hp.run()