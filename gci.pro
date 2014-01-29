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
    gciString.cpp \
    gciWavefunction.cpp

HEADERS += \
    gciHamiltonian.h \
    gciParameters.h \
    gciNode.h \
    gciExcitation.h \
    gciDeterminant.h \
    gciString.h \
    gciWavefunction.h

html.target = $$PWD/html
html.commands = (cd $$PWD ; doxygen $$PWD/Doxyfile;)
html.depends = $(OBJECTS)
dox.target = html
dox.depends = $$PWD/html
QMAKE_EXTRA_TARGETS += dox html all

all.target = all
all.depends += FCIDUMP
all.depends += html

input.target = FCIDUMP
input.commands = ln -s $$PWD/FCIDUMP FCIDUMP
input.depends =
QMAKE_EXTRA_TARGETS += input






