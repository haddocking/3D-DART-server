# 3D-DART-server

3D-DART web server repository. Docker is needed in order to build the 3D-DART server container and to execute it. It is based in the Ubuntu 18.04 Docker image (with 32bit support needed for the version of X3DNA software used by 3D-DART) and with Apache2 and Python 2.7.

## Build Docker container

```bash
cd 3D-DART-server
docker build -t 3d-dart .
```

## Run 3D-DART Docker container

```bash
docker run --name 3ddart -p 80:80 -i -t 3d-dart
```


