n_threads=4
local_dir=/mnt/research/ged/camille/lamprey/2013-lamprey-data/

environment:
	conda env export --file environment.conda

list:
	./pipeline list --resources-metadata resources.json --config-metadata config.json --local-file-dir $(local_dir)

busco:
	./pipeline busco --resources-metadata resources.json --config-metadata config.json --local-file-dir $(local_dir)

databases:
	./pipeline databases --resources-metadata resources.json --config-metadata config.json --local-file-dir $(local_dir)

clean:
	./pipeline clean --resources-metadata resources.json --config-metadata config.json --local-file-dir $(local_dir)
