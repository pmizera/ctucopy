# You should modify the man page and then create html, pdf files.
man2html  ./ctucopy4 > ctucopy_400.html
man -t ./ctucopy4 |  ps2pdf -  ctucopy_400.pdf
