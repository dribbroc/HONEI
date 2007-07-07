CLEANFILES = *~ cleanauthors cleancontributors

MAINTAINERCLEANFILES = Makefile.in

doxygen: doxygen.conf htmlpages
	doxygen $(top_srcdir)/doc/doxygen.conf
	cp powered-by-*.png $(top_srcdir)/doc/html/

HTMLPAGES = \
	authors.html

htmlpages: $(HTMLPAGES)
	for f in $(HTMLPAGES) pg512.css ; do cp $(srcdir)/$$f $(srcdir)/html/ ; done

cleanauthors: $(top_srcdir)/AUTHORS
	sed -e '1,3d' \
	    -e '/^Contributions/,$$d' \
	    -e 's,&,\&amp;,g' -e 's,<,\&lt;,g' -e 's,>,\&gt;,g' -e 's,@, at ,g' \
	    -e 's,^    \(.*\),<dd>\1</dd>,' \
	    -e 's,^\([^ <].*\),<dt>\1</dt>,' \
	    < ../AUTHORS > cleanauthors

cleancontributors: $(top_srcdir)/AUTHORS
	sed -e '1,/^Contributions/d' \
	    -e '/^=*$$/d' \
	    -e 's,&,\&amp;,g' -e 's,<,\&lt;,g' -e 's,>,\&gt;,g' -e 's,@, at ,g' \
	    -e 's,^    \(.*\),<dd>\1</dd>,' \
	    -e 's,^\([^ <].*\),<dt>\1</dt>,' \
	    < ../AUTHORS > cleancontributors

authors.html: authors.html.skel cleanauthors cleancontributors
	sed -e '/@AUTHORS@/r cleanauthors' -e '/@AUTHORS@/d' \
	    -e '/@CONTRIBUTORS@/r cleancontributors' -e '/@CONTRIBUTORS@/d' \
	    -e '/@HEADER@/r header.html' \
	    -e '/@HEADER@/d' \
	    -e '/@FOOTER@/r footer.html' \
	    -e '/@FOOTER@/d' < $(srcdir)/authors.html.skel > authors.html



