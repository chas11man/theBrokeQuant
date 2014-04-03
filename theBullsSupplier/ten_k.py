import json
import urllib2


with open('dumps/api.json', 'w') as f:
	json.dump(json.load(urllib2.urlopen('http://api.corpwatch.org/2008/companies.json?company_name=google')), f)
