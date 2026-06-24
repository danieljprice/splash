FROM debian:12.11

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    gfortran \
    git \
    make \
    cmake \
    libpng-dev \
    libx11-dev \
    libcairo2-dev \
    pkg-config \
    vim \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/splash
RUN git clone https://github.com/danieljprice/splash.git . && \
    git clone https://github.com/danieljprice/giza.git giza && \
    git checkout v3.12.0 && \
    make SYSTEM=gfortran withgiza

ENV PATH="/opt/splash/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/splash/giza/lib:$LD_LIBRARY_PATH"

RUN pip install pillow "imagehash>=4.3.2" "numpy<2.0" \
    --break-system-packages

WORKDIR /workspace