## Process this file with automake to generate Makefile.in

SUBDIRS = src include

export: dist
	scp $(PACKAGE)-$(VERSION).tar.gz phsaap@sirius-b.csc.warwick.ac.uk:/data7/phsaap/software/.
	ssh phsaap@sirius-b.csc.warwick.ac.uk "cd /data7/phsaap/software ; \rm -f $(PACKAGE).tar.gz ; ln -s $(PACKAGE)-$(VERSION).tar.gz $(PACKAGE).tar.gz; \rm -fr $(PACKAGE); mkdir $(PACKAGE)"

