VERSION_MAJOR = 1
VERSION = 1.0.0
PROJECT_NAME=biquad
PREFIX=/usr/local

LIB_NAME = lib$(PROJECT_NAME).so
SONAME = lib$(PROJECT_NAME).so.$(VERSION_MAJOR)

HEADERS = biquad.h

OBJECTS = biquad.o

TARGET_ARCH := $(shell uname -m)

CFLAGS = -fPIC -O2

LIBS = -lm

TARGETS = $(SONAME) $(LIB_NAME)

$(SONAME): $(LIB_NAME)
	ln -sf $(LIB_NAME) $(SONAME)

all: $(TARGETS)

install: $(LIB_NAME)
	cp $(LIB_NAME) $(PREFIX)/lib/$(LIB_NAME).$(VERSION)
	ln -sf $(LIB_NAME).$(VERSION) $(PREFIX)/lib/$(LIB_NAME)
	ldconfig

$(LIB_NAME): $(OBJECTS)
	$(CC) -shared -Wl,-soname,$(SONAME) -o $@ $(OBJECTS) $(LIBS)

%.o: %.c
	$(CC) -c -o $@ $(CFLAGS) $(LIBS) $<

test: test.c $(OBJECTS)
	$(CC) -o $@ $(CFLAGS) test.c -l$(PROJECT_NAME) $(LIBS)

clean:
	rm -f $(TARGETS) $(OBJECTS)

