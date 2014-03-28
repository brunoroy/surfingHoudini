QMAKE_CXXFLAGS += -std=gnu++0x

HOUDINI_VERSION = "12.5.562"

TARGET = SurfingSAT
TEMPLATE = lib
CONFIG += plugin no_plugin_name_prefix

HEADERS += \
    include/*.h

SOURCES += \
    build/sesitag.C \
    src/*.cpp

DEFINES += \
   USE_PTHREADS  \
   GCC3 \
   SIZEOF_VOID_P=8 \
   DLLEXPORT="" \
   _GNU_SOURCE \
   ENABLE_THREADS \
   ENABLE_UI_THREADS \
   VERSION=HOUDINI_VERSION \
   LINUX \
   AMD64 \
   MAKING_DSO \
   SESI_LITTLE_ENDIAN \
   SWAP_BITFIELDS

INCLUDEPATH += include external /opt/hfs$$HOUDINI_VERSION/toolkit/include /opt/hfs$$HOUDINI_VERSION/toolkit/include/htools /opt/hfs$$HOUDINI_VERSION/toolkit/include/OpenEXR

#Houdini toolkit
LIBS += -L/opt/hfs$$HOUDINI_VERSION/dsolib -lHoudiniUI -lHoudiniOPZ \
        -lHoudiniOP3 -lHoudiniOP2 -lHoudiniOP1 \
        -lHoudiniSIM -lHoudiniGEO -lHoudiniPRM -lHoudiniUT \
        -ldwelf

#libboost
LIBS += -L/opt/hfs$$HOUDINI_VERSION/dsolib -lboost_iostreams -lboost_regex -lboost_thread -lboost_date_time -lboost_system

#libtbb
LIBS += -L/opt/hfs$$HOUDINI_VERSION/dsolib -ltbb

#OpenCL
LIBS += -L/opt/hfs$$HOUDINI_VERSION/dsolib -lOpenCL

#others
LIBS += -lImath

CONFIG(debug, debug|release) {
    DESTDIR = build
    OBJECTS_DIR = build
    MOC_DIR = build
    UI_DIR = build
} else {
    DESTDIR = build
    OBJECTS_DIR = build
    MOC_DIR = build
    UI_DIR = build
}
