#-------------------------------------------------
#
# Project created by QtCreator 2012-03-06T06:11:54
#
#-------------------------------------------------

QT       -= core

QT       -= gui

TARGET = treeci
CONFIG   += console
CONFIG   -= app_bundle

#INCLUDEPATH += /usr/local/boost

TEMPLATE = app


SOURCES += main.cpp \
    TreeCIHamiltonian.cpp \
    TreeCIParameters.cpp \
    TreeCINode.cpp \
    TreeCIExcitation.cpp \
    TreeCIDeterminant.cpp \
    TreeCIString.cpp

HEADERS += \
    TreeCIHamiltonian.h \
    TreeCIParameters.h \
    TreeCINode.h \
    TreeCIExcitation.h \
    TreeCIDeterminant.h \
    TreeCIString.h





