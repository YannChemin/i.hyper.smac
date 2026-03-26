GISBASE ?= $(shell grass --config path)

CC      = gcc
CFLAGS  = -O2 -fPIC -Wall \
          -I$(GISBASE)/include
LDFLAGS = -shared -L$(GISBASE)/lib \
          -lgrass_g3d -lgrass_gis \
          -Wl,-rpath,$(GISBASE)/lib

TARGET = libr3dslice.so
SRC    = r3d_slice.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o

.PHONY: all clean
