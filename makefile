all: create grasper

clean:
	rm src/classes/*.class
	rmdir src/classes
	rm -rf bin/*.jar
	rmdir bin
create:
	mkdir -p bin

grasper:
#	javac -Xlint:unchecked Thread.java 
	mkdir -p src/classes
	javac -cp src src/Grasper.java -d src/classes
	jar cfe bin/grasper.jar Grasper -C src/classes .
