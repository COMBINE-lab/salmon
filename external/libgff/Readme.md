LibGFF
======

This is an attempt to perform a simple "libraryfication" of the GFF/GTF parsing
code that is used in the [Cufflinks](http://cufflinks.cbcb.umd.edu/index.html)
codebase.  There are not many (any?) relatively lightweight GTF/GFF parsers
exposing a C++ interface, and the goal of this library is to provide this
functionality without the necessity of drawing in a heavy-weight dependency
like SeqAn.  
