CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		Alignment.o	BandedGlobalAlignment.o

LIBS =

TARGET =	kgalign.exe

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
