#----------------------------------------------------------------
#   Parent Makefile for SPLASH
#   This file is just a wrapper for the sub-make in the build
#   directory. Refer to build/Makefile for more details.
#
#   (c) 2007-2013 Daniel Price
#
#----------------------------------------------------------------

.PHONY: splash install docs tests src bin
splash:
	@cd build; ${MAKE} ${MAKECMDGOALS}

%::
	@cd build; ${MAKE} ${MAKECMDGOALS}

install:
	@cd build; ${MAKE} ${MAKECMDGOALS}

docs:
	@cd build; ${MAKE} ${MAKECMDGOALS}

clean:
	@cd build; ${MAKE} clean
