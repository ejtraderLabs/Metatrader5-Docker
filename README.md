## Usage

Metatrader 5 on docker and VNC
The container runs a VNC server on port 5900. This port has to be mapped for VNC clients to access it:
cd
vnc login: root password: root

run and build image named as ejtraderMT and run container named as ejtraderMT

```bash
make run
```

build image named ejtraderMT

```bash
make build
```

login to shell

```bash
 make shell
```
