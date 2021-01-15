
include Makefile.in

#############################################################################

INCLUDE = -I $(LIB_DIR)
RPATH = -rpath=$(LIB_DIR)

COPTS += $(INCLUDE)
CXXOPTS += $(INCLUDE)

###### targets ##############################################################

### (the main target makes use of the "libINtfit.so" library)
all: lib
	$(CXX) $(CXXOPTS) $(DEBUG) -I $(LIB_DIR) main.cpp \
         -Wl,$(RPATH) -L $(LIB_DIR) -lINtfit -lm

### (the software is made into a library: "libINtfit.so"
lib:
	$(CXX) $(CXXOPTS) $(DEBUG) -I $(LIB_DIR) -c intfit.cpp
	$(CXX) -shared -Wl,-rpath=$(LIB_DIR),-soname,libINtfit.so -o libINtfit.so \
          intfit.o -lc -lm 

clean:
	rm -f *.o a.out *.so 

