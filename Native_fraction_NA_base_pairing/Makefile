#
# Usage:
# $ source /path/to/GMXRC
# $ make -f Makefile.pkg
#
#change the name of the program here
NAME=g_native_hb

#add extra c file to compile here
EXTRA_SRC=

###############################################################3
#below only boring default stuff
#only change it if you know what you are doing ;-)

#what should be done by default
all: $(NAME)

#if GMXLDLIB is defined we add it to PKG_CONFIG_PATH
ifeq "$(origin GMXLDLIB)" "undefined"
  $(error "GMXLDLIB not found, please source GMXRC")
else
  export PKG_CONFIG_PATH:=${PKG_CONFIG_PATH}:${GMXLDLIB}/pkgconfig
  export PKG_CONFIG_PATH:=${PKG_CONFIG_PATH}:/Users/guillem/work/soft/fftw3/lib/pkgconfig
endif

#get CPPFLAGS and LDFLAGS from pkg-config
CPPFLAGS=`pkg-config --cflags libgmx libmd ` -O3
LDFLAGS=`pkg-config --libs libgmx  libmd ` 

#generate a list of object (.o) files
OBJS=$(patsubst %.c,%.o,$(NAME).c $(EXTRA_SRC))

#main program depend on all objects, rest is done by implicit rules
$(NAME): $(OBJS)

#clean up rule
clean:
	rm -f $(NAME) $(OBJS)

#all, clean are phony rules, e.g. they are always run
.PHONY: all clean
