#Stole template from here:
#http://www.itnotes.de/docker/development/tools/2014/08/31/speed-up-your-docker-workflow-with-a-makefile/
NS = humanlongevity
TAG = latest

BASENAME = splithunter

.PHONY: all
all: build push
test: build testrun

build:
	docker system prune -f
	docker build -f $(BASENAME).dockerfile -t $(NS)/$(BASENAME):$(TAG) ./
push:
	docker push $(NS)/$(BASENAME):$(TAG)
testrun:
	docker run $(NS)/$(BASENAME):$(TAG) splithunter.py
