## Process this file with automake to generate Makefile.in

SUBDIRS = src include

export: dist
	scp $(PACKAGE)-$(VERSION).tar.gz $(WEB_SERVE):$(WEB_PATH)/software/.
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software ; \rm -f $(PACKAGE).tar.gz ; ln -s $(PACKAGE)-$(VERSION).tar.gz $(PACKAGE).tar.gz; \rm -fr $(PACKAGE); mkdir $(PACKAGE)"
	scp -r doc/html $(WEB_SERVE):$(WEB_PATH)/software/$(PACKAGE)/.
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software; sed -e 's%<\!-- mem version --><td>.*</td>%<\!-- mem version --><td>$(VERSION)</td>%' index.html > junk.html; mv junk.html index.html"
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software; sed -e 's%<\!-- mem date --><td>.*</td>%<\!-- mem date --><td>$(DATE)</td>%' index.html > junk.html; mv junk.html index.html"


