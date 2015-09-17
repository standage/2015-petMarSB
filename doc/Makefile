BIBFILE=petmarrna.bib

all: petmarrna.draft.pdf petmarrna.draft.rst

petmarrna.draft.rst: petmarrna.draft.tex
	pandoc $< -o $@

petmarrna.draft.pdf: petmarrna.draft.tex
	rm -rf tex_out
	mkdir -p tex_out
	cd tex_out; ln -s ../$<
	cd tex_out; ln -s ../$(BIBFILE)
	cd tex_out; pdflatex `basename -s .tex $<`; bibtex `basename -s .tex $<`; pdflatex `basename -s .tex $<`; pdflatex `basename -s .tex $<`
	cp tex_out/$@ .

clean:
	rm -rf tex_out petmarrna.draft.pdf petmarrna.draft.tex
