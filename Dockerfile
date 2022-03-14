FROM python:latest

RUN pip install haploblock-shuffler==0.0.4
RUN apt update && apt install -y tabix
