all: create samsort clusterapp

clean:
	rm *.class
	rm -rf bin/*.jar

create:
	mkdir -p bin

samsort:
	javac SAMSort.java
	jar cfe bin/SAMSort.jar SAMSort -C . .

clusterapp:
	javac ClusterApp.java
	jar cfe bin/ClusterApp.jar ClusterApp -C . .
