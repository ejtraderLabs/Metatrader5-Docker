build:
	docker build -t feed .

run: build
	docker run --rm -dit -p 5900:5900 -p 15555:15555 -p 15556:15556 -p 15557:15557 -p 15558:15558 --name ejtraderMT -v ejtraderMT:/data ejtraderMT

shell: 
	docker exec -it datafeed bash

users: build
	docker exec -it datafeed adduser novouser