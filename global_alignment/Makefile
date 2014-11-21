CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		Alignment.o	GlobalAlignment.o

LIBS =

TARGET =	galign.exe

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
