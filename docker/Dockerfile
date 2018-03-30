#---------------------------
FROM rootproject/root-ubuntu16:6.10

USER root

RUN apt-get update && \
  apt-get -y install \
  vim \
  nano \
  libgsl2 \
  python-numpy \
  perl \
  bc \
  git \
  liblapack3

RUN ln -s /usr/lib/x86_64-linux-gnu/libgsl.so.19.0.0 /usr/lib/libgsl.so.0 && \
  echo '/usr/local/bin/root -b $@' > /usr/local/sbin/root && chmod +x /usr/local/sbin/root &&\
  git clone https://github.com/acidbourbon/mdc_3d_garfield.git

#---------------------------