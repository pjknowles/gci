#-------------------------------------------------
#
# Project created by QtCreator 2012-03-06T06:11:54
#
#-------------------------------------------------

QT       -= core

QT       -= gui

TARGET = gci.exe
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    gciHamiltonian.cpp \
    gciParameters.cpp \
    gciNode.cpp \
    gciExcitation.cpp \
    gciDeterminant.cpp \
    gciString.cpp

HEADERS += \
    gciHamiltonian.h \
    gciParameters.h \
    gciNode.h \
    gciExcitation.h \
    gciDeterminant.h \
    gciString.h

dox.target = doxygen
dox.commands = doxygen $$PWD/Doxyfile;
dox.depends =
QMAKE_EXTRA_TARGETS += dox

input.target = input
input.commands = ln -s $$PWD/FCIDUMP FCIDUMP
input.depends =
QMAKE_EXTRA_TARGETS += input






