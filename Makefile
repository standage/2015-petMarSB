n_threads=4

environment:
	conda env export --file environment.conda

list:
	./pipeline list --resources-metadata resources.json --config-metadata config.json

print-tasks:
	./pipeline --print-tasks --resources-metadata resources.json --config-metadata config.json

busco:
	./pipeline busco --resources-metadata resources.json --config-metadata config.json

blast:
	./pipeline blast --resources-metadata resources.json --config-metadata config.json

databases:
	./pipeline databases --resources-metadata resources.json --config-metadata config.json

test:
	./pipeline run --resources-metadata test/resources.json --config-metadata test/config.json --local-file-dir test/ --assembly-file lamp-test.fasta

clean:
	./pipeline clean --resources-metadata resources.json --config-metadata config.json
