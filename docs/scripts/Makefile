.PHONY: default clean


default: doxyCFG_TCSAM02
	./tpl2dox.sh ../../tcsam02.tpl
	doxygen doxyCFG_TCSAM02


clean:
	@cd html; rm -rf d* *.* search
