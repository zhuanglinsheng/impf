FROM python:3.9-slim

RUN apt-get update && \
    apt-get install -y valgrind build-essential

# RUN apt-get install -y <Libraries>

WORKDIR /app

COPY . /app

CMD ["bash"]