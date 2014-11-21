CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		Alignment.o	MultipleGlobalAlignment.o

LIBS =

TARGET =	malign.exe

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
