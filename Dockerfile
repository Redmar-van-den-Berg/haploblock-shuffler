FROM python:latest

RUN apt update && apt install -y tabix
RUN pip install haploblock-shuffler==0.0.6
