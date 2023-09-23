FROM ubuntu:22.04

RUN apt-get update
RUN apt-get install -y curl

RUN curl -L -O https://github.com/marschall-lab/GFAffix/releases/download/0.1.5/GFAffix-0.1.5_linux_x86_64.tar.gz

RUN tar -xzvf GFAffix-0.1.5_linux_x86_64.tar.gz

RUN mv ./GFAffix-0.1.5_linux_x86_64/gfaffix /usr/local/bin/
