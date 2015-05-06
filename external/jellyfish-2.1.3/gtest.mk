##############################
# Gtest build.
##############################
# Build rules for libraries.
check_LTLIBRARIES = libgtest.la libgtest_main.la

libgtest_la_SOURCES = unit_tests/gtest/src/gtest-all.cc
libgtest_main_la_SOURCES = unit_tests/gtest/src/gtest_main.cc
libgtest_main_la_LIBADD = libgtest.la
libgtest_la_CXXFLAGS = -I$(top_srcdir)/unit_tests
libgtest_main_la_CXXFLAGS = -I$(top_srcdir)/unit_tests

GTEST_SRC = unit_tests/gtest/src/gtest-all.cc	\
	    unit_tests/gtest/src/gtest_main.cc	\
	    unit_tests/gtest/gtest.h

EXTRA_DIST += $(GTEST_SRC)

