CONFIG        += c++11
QMAKE_CXXFLAGS += -Werror
HEADERS       = window.h \
                parab.h \
                bessel.h
SOURCES       = main.cpp \
                window.cpp \
                parab.cpp \
                bessel.cpp
QT += widgets
