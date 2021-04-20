build:
	docker build -t ejtraderMT .

run: build
	docker run --rm -dit -p 5900:5900 -p 15555:15555 -p 15556:15556 -p 15557:15557 -p 15558:15558 --name ejtraderMT -v ejtraderMT:/data ejtraderMT

shell: 
	docker exec -it ejtraderMT bash

users: build
	docker exec -it ejtraderMT adduser novouser