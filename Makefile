default: all

all: ${FILES}
	$(MAKE) -C lib build

clean:
	$(MAKE) -C lib clean

