# docker build --build-arg ver=0.1.0

FROM python:3.8-slim

RUN apt-get update -y -m && \
    DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
                                                   gcc \
                                                   g++ && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN python -m pip install --upgrade pip && \
    python -m pip install deltapd==0.1.1

ENTRYPOINT ["deltapd"]
CMD ["--help"]
