n_threads=4

all: doc/petmarrna.draft.pdf

doc/petmarrna.draft.tex: petmarrna.draft.texw
	Pweave -f tex $< -s ipython; mv petmarrna.draft.tex doc/

doc/petmarrna.draft.pdf: doc/petmarrna.draft.tex
	cd doc; make

environment: FORCE
	conda env export --file environment.conda

list: FORCE
	./pipeline list --resources-metadata resources.json --config-metadata config.json

print-tasks: FORCE
	./pipeline --print-tasks --resources-metadata resources.json --config-metadata config.json

busco: FORCE
	./pipeline busco --resources-metadata resources.json --config-metadata config.json

blast: FORCE
	./pipeline blast --resources-metadata resources.json --config-metadata config.json

databases: FORCE
	./pipeline databases --resources-metadata resources.json --config-metadata config.json

test: FORCE
	./pipeline run --resources-metadata test/resources.json --config-metadata test/config.json --local-file-dir test/ --assembly-file lamp-test.fasta

clean:
	cd doc; make clean

FORCE:
