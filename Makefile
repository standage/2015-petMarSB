n_threads=4

environment:
	conda env export --file environment.conda

list:
	./pipeline list --resources-metadata resources.json --config-metadata config.json

busco:
	./pipeline busco --resources-metadata resources.json --config-metadata config.json


clean:
	./pipeline clean --resources-metadata resources.json --config-metadata config.json
